#include "UncollidedFluxRayStudy.h"

#include "RealSphericalHarmonics.h"

#include "libmesh/parallel_algebra.h"

registerMooseObject("GnatApp", UncollidedFluxRayStudy);

// #define DEBUG_OUTPUT

InputParameters
UncollidedFluxRayStudy::validParams()
{
  auto params = RayTracingStudy::validParams();
  params.addClassDescription(
      "A ray tracing study which generates rays from point, surface, and volume sources for "
      "computing the uncollided component of the angular flux.");

  params.addRequiredParam<unsigned int>(
      "num_groups",
      "The number of spectral energy groups that this study is computing the uncollided flux for.");

  // Point sources.
  params.addParam<std::vector<Point>>("point_source_locations",
                                      "The locations of all point sources in the problem space.");
  params.addParam<std::vector<std::vector<Real>>>(
      "point_source_moments",
      "A double vector containing a list of external source moments for all point particle "
      "sources. The external vector should correspond with the order of "
      "'point_source_locations'.");
  params.addParam<std::vector<unsigned int>>(
      "point_source_anisotropies",
      "The anisotropies of the point sources. The vector should correspond with the order of "
      "'point_source_locations'");

  // Surface sources.
  params.addParam<std::vector<BoundaryName>>("source_boundaries",
                                             "The boundaries to apply incoming "
                                             "flux boundary conditions.");
  params.addParam<std::vector<std::vector<Real>>>(
      "boundary_source_moments",
      "A double vector containing the external source moments for "
      "all boundaries. The exterior vector must correspond with the surface source boundary "
      "conditions provided in 'source_boundaries'.");
  params.addParam<std::vector<unsigned int>>(
      "boundary_source_anisotropy",
      "The degree of anisotropy of the boundary source moments. The exterior vector must "
      "correspond with the surface source boundary "
      "conditions provided in 'source_boundaries'.");

  // Volume sources.
  params.addParam<std::vector<SubdomainName>>("volumetric_source_blocks",
                                              "The list of blocks (ids or "
                                              "names) that host a volumetric source.");
  params.addParam<std::vector<std::vector<Real>>>(
      "volumetric_source_moments",
      "A double vector containing a list of external source moments for all volumetric particle "
      "sources. The external vector should correspond with the order of "
      "'volumetric_source_blocks'.");
  params.addParam<std::vector<unsigned int>>(
      "volumetric_source_anisotropies",
      "The anisotropies of the volumetric sources. The vector should correspond with the order of "
      "'volumetric_source_blocks'");

  // Spatial quadrature parameters.
  MooseEnum qorders("CONSTANT FIRST SECOND THIRD FOURTH FIFTH SIXTH SEVENTH EIGHTH NINTH TENTH "
                    "ELEVENTH TWELFTH THIRTEENTH FOURTEENTH FIFTEENTH SIXTEENTH SEVENTEENTH "
                    "EIGHTTEENTH NINTEENTH TWENTIETH",
                    "CONSTANT");
  params.addParam<MooseEnum>("volume_order",
                             qorders,
                             "The volume quadrature rule order. For simplicity the same quadrature "
                             "order is used for volumetric sources and target elements.");
  params.addParam<MooseEnum>("face_order", qorders, "The face quadrature rule order.");

  MooseEnum qtypes("GAUSS GRID", "GAUSS");
  params.addParam<MooseEnum>("volume_type",
                             qtypes,
                             "The volume quadrature type. For simplicity the same quadrature type "
                             "is used for both volumetric sources and target elements.");
  params.addParam<MooseEnum>("face_type", qtypes, "The face quadrature type.");

  params.addParam<std::string>("source_and_weights_name",
                               "source_and_weights",
                               "The name of the ray data which houses the source intensity "
                               "and spatial weights, pre-multiplied.");

  //----------------------------------------------------------------------------
  // Angular quadrature parameters.
  params.addRangeCheckedParam<unsigned int>("n_polar",
                                            30,
                                            "n_polar > 0",
                                            "Number of Legendre polar "
                                            "quadrature points in a single "
                                            "octant of the unit sphere. "
                                            "Defaults to 3.");
  params.addRangeCheckedParam<unsigned int>("n_azimuthal",
                                            30,
                                            "n_azimuthal > 0",
                                            "Number of Chebyshev azimuthal "
                                            "quadrature points in a single "
                                            "octant of the unit sphere. "
                                            "Defaults to 3.");

  // It's impractical to register rays due to the sheer number of them.
  params.set<bool>("_use_ray_registration") = false;

  // Run this object when a timestep begins.
  params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_BEGIN;

  return params;
}

UncollidedFluxRayStudy::UncollidedFluxRayStudy(const InputParameters & parameters)
  : RayTracingStudy(parameters),
    _target_in_element(registerRayData("target_source_same_element")),
    _dim(_mesh.dimension()),
    _num_groups(getParam<unsigned int>("num_groups")),
    _symmetry_factor(_dim == 2u ? 2.0 : 1.0),
    _volume_fe(FEBase::build(_dim, FEType(CONSTANT, MONOMIAL))),
    _q_volume(QBase::build(Moose::stringToEnum<QuadratureType>(getParam<MooseEnum>("volume_type")),
                           _dim,
                           Moose::stringToEnum<Order>(getParam<MooseEnum>("volume_order")))),
    _face_fe(FEBase::build(_dim, FEType(CONSTANT, MONOMIAL))),
    _q_face(QBase::build(Moose::stringToEnum<QuadratureType>(getParam<MooseEnum>("face_type")),
                         _dim - 1u,
                         Moose::stringToEnum<Order>(getParam<MooseEnum>("face_order")))),
    _2D_angular_quadrature(
        _dim == 2u ? std::make_unique<LegendrePolynomial>(4u * getParam<unsigned int>("n_polar"))
                   : nullptr),
    _3D_angular_quadrature(_dim == 3u ? std::make_unique<GaussAngularQuadrature>(
                                            2u * getParam<unsigned int>("n_azimuthal"),
                                            2u * getParam<unsigned int>("n_polar"),
                                            MajorAxis::X,
                                            ProblemType::Cartesian3D)
                                      : nullptr),
    _num_dir(0u),
    _point_source_locations(getParam<std::vector<Point>>("point_source_locations")),
    _point_source_moments(getParam<std::vector<std::vector<Real>>>("point_source_moments")),
    _point_source_anisotropy(getParam<std::vector<unsigned int>>("point_source_anisotropies")),
    _source_boundary_names(getParam<std::vector<BoundaryName>>("source_boundaries")),
    _boundary_source_moments(getParam<std::vector<std::vector<Real>>>("boundary_source_moments")),
    _boundary_source_anisotropy(getParam<std::vector<unsigned int>>("boundary_source_anisotropy")),
    _volume_source_blocks(getParam<std::vector<SubdomainName>>("volumetric_source_blocks")),
    _volume_source_moments(getParam<std::vector<std::vector<Real>>>("volumetric_source_moments")),
    _volume_source_anisotropy(getParam<std::vector<unsigned int>>("volumetric_source_anisotropies"))

{
  _volume_fe->attach_quadrature_rule(_q_volume.get());
  _face_fe->attach_quadrature_rule(_q_face.get());

  _volume_fe->get_xyz();
  _face_fe->get_xyz();

  _source_spatial_weights.reserve(_num_groups);
  for (unsigned int g = 0u; g < _num_groups; ++g)
    _source_spatial_weights.emplace_back(registerRayData(
        getParam<std::string>("source_and_weights_name") + "_" + Moose::stringify(g)));

  if (_dim <= 1u)
    mooseError("Ray tracing for the uncollided flux is not supported on 1D meshes as ray effects "
               "don't exist in 1D.");

  // Handle possible errors for point sources.
  if (_point_source_locations.size() > 0u)
  {
    if (_point_source_locations.size() != _point_source_moments.size() ||
        _point_source_locations.size() != _point_source_anisotropy.size())
      mooseError(
          "Mismatch between the number of point sources and the provided anisotropy / moments.");

    for (unsigned int i = 0u; i < _point_source_locations.size(); ++i)
    {
      const auto num_moments =
          _dim == 3u ? (_point_source_anisotropy[i] + 1u) * (_point_source_anisotropy[i] + 1u)
                     : (_point_source_anisotropy[i] + 1u) * (_point_source_anisotropy[i] + 2u) / 2u;

      if (num_moments * _num_groups != _point_source_moments[i].size())
        mooseError("Mismatch between the number of provided moments for the point source at " +
                   Moose::stringify(i) + " with the number of provided groups / anisotropy.");
    }
  }

  // Handle possible errors for surface sources.
  if (_source_boundary_names.size() > 0u)
  {
    if (_source_boundary_names.size() != _boundary_source_moments.size() ||
        _source_boundary_names.size() != _boundary_source_anisotropy.size())
      mooseError("Mismatch between the number of surface sources and the provided anisotropy / "
                 "moments.");

    for (unsigned int i = 0u; i < _source_boundary_names.size(); ++i)
    {
      const auto num_moments =
          _dim == 3u
              ? (_boundary_source_anisotropy[i] + 1u) * (_boundary_source_anisotropy[i] + 1u)
              : (_boundary_source_anisotropy[i] + 1u) * (_boundary_source_anisotropy[i] + 2u) / 2u;

      if (num_moments * _num_groups != _boundary_source_moments[i].size())
        mooseError("Mismatch between the number of provided moments for the surface source at " +
                   Moose::stringify(i) + " with the number of provided groups / anisotropy.");
    }
  }

  // Handle possible errors for volumetric sources.
  if (_volume_source_blocks.size() > 0u)
  {
    if (_volume_source_blocks.size() != _volume_source_moments.size() ||
        _volume_source_blocks.size() != _volume_source_anisotropy.size())
      mooseError("Mismatch between the number of volumetric sources and the provided anisotropy / "
                 "moments.");

    for (unsigned int i = 0u; i < _volume_source_blocks.size(); ++i)
    {
      const auto num_moments =
          _dim == 3u
              ? (_volume_source_anisotropy[i] + 1u) * (_volume_source_anisotropy[i] + 1u)
              : (_volume_source_anisotropy[i] + 1u) * (_volume_source_anisotropy[i] + 2u) / 2u;

      if (num_moments * _num_groups != _volume_source_moments[i].size())
        mooseError("Mismatch between the number of provided moments for the volumetric source at " +
                   Moose::stringify(i) + " with the number of provided groups / anisotropy.");
    }
  }

  _num_dir = _dim == 2u ? _2D_angular_quadrature->degree() : _3D_angular_quadrature->totalOrder();
}

void
UncollidedFluxRayStudy::cartesianToSpherical(const RealVectorValue & direction,
                                             Real & mu,
                                             Real & omega)
{
  mu = 0.0;
  omega = 0.0;

  mu = direction(0);

  if (direction(1) > 0.0)
  {
    if (direction(2) > 0.0)
      omega += std::atan(std::abs(direction(2)) / std::abs(direction(1)));
    else if (direction(2) < 0.0)
      omega += std::atan(std::abs(direction(2)) / std::abs(direction(1))) + libMesh::pi;
  }
  else if (direction(1) < 0.0)
  {
    if (direction(2) > 0.0)
      omega += std::atan(std::abs(direction(2)) / std::abs(direction(1))) + 3.0 * libMesh::pi / 2.0;
    else if (direction(2) < 0.0)
      omega += std::atan(std::abs(direction(2)) / std::abs(direction(1))) + libMesh::pi / 2.0;
    else
      omega += libMesh::pi;
  }
  else
  {
    if (direction(2) > 0.0)
      omega += libMesh::pi / 2.0;
    else if (direction(2) < 0.0)
      omega += 3.0 * libMesh::pi / 2.0;
  }
}

Real
UncollidedFluxRayStudy::computeSHSource(unsigned int source,
                                        const std::vector<std::vector<Real>> & moments,
                                        const std::vector<unsigned int> & anisotropy,
                                        const RealVectorValue & direction,
                                        unsigned int group_index)
{
  const auto num_moments = _dim == 3u ? (anisotropy[source] + 1u) * (anisotropy[source] + 1u)
                                      : (anisotropy[source] + 1u) * (anisotropy[source] + 2u) / 2u;

  unsigned int moment_index = group_index * num_moments;

  Real mu = 0.0;
  Real omega = 0.0;
  cartesianToSpherical(direction, mu, omega);

  Real src_l = 0.0;
  Real src_m = 0.0;
  if (_dim == 2u)
  {
    for (unsigned int l = 0u; l <= anisotropy[source]; ++l)
    {
      for (int m = 0; m <= static_cast<int>(l); ++m)
      {
        src_m += moments[source][moment_index] * RealSphericalHarmonics::evaluate(l, m, mu, omega);
        moment_index++;
      }

      src_l += src_m * (2.0 * static_cast<Real>(l) + 1.0) / (4.0 * libMesh::pi) * _symmetry_factor;
      src_m = 0.0;
    }
  }
  else
  {
    for (unsigned int l = 0u; l <= anisotropy[source]; ++l)
    {
      for (int m = -1 * static_cast<int>(l); m <= static_cast<int>(l); ++m)
      {
        src_m += moments[source][moment_index] * RealSphericalHarmonics::evaluate(l, m, mu, omega);
        moment_index++;
      }

      src_l += src_m * (2.0 * static_cast<Real>(l) + 1.0) / (4.0 * libMesh::pi) * _symmetry_factor;
      src_m = 0.0;
    }
  }

  return src_l;
}

void
UncollidedFluxRayStudy::generateRays()
{
  // Tally the total number of rays generated along with the number of source and target points.
  std::size_t total_num_rays = 0u;
  std::size_t num_point_source_points = 0u;
  std::size_t num_surface_source_points = 0u;
  std::size_t num_volume_source_points = 0u;

#ifdef DEBUG_OUTPUT
  std::size_t debug_num_rays = 0u;
#endif

  // Generate the global quadrature points.
  std::vector<Point> global_spatial_q_points;
  std::vector<Real> global_spatial_q_weights;
  {
    const auto & points = _volume_fe->get_xyz();
    const auto & weights = _volume_fe->get_JxW();

    // Gather a list of global quadrature points and weights.
    for (const auto & elem : *_mesh.getActiveLocalElementRange())
    {
      _volume_fe->reinit(elem);
      for (unsigned int i = 0u; i < points.size(); ++i)
      {
        global_spatial_q_points.emplace_back(points[i]);
        global_spatial_q_weights.emplace_back(weights[i]);
      }
    }

    // Send all of the points/weights for the targets to each processor.
    _comm.allgather(global_spatial_q_points);
    _comm.allgather(global_spatial_q_weights);
  }

  // Point sources.
  if (_point_source_locations.size() > 0u)
  {
    // Loop over point sources to find the element they live in.
    std::unordered_map<unsigned int, const Elem *> local_point_elements;

    const auto locator = _mesh.getPointLocator();
    for (unsigned int i = 0u; i < _point_source_locations.size(); ++i)
    {
      const auto elem = (*locator)(_point_source_locations[i]);
      if (!elem)
        mooseWarning("The point source with index " + Moose::stringify(i) +
                     " does not exist on the mesh!");

      // This processor owns this starting element.
      if (elem->processor_id() == _pid)
        local_point_elements.emplace(i, elem);
    }

    total_num_rays += local_point_elements.size() * global_spatial_q_points.size();
    num_point_source_points += local_point_elements.size();
    reserveRayBuffer(local_point_elements.size() * global_spatial_q_points.size());

    // Now that we've found the point sources we own, we set up the rays.
    for (unsigned int i = 0u; i < global_spatial_q_points.size(); ++i)
    {
      for (const auto & [src_index, elem] : local_point_elements)
      {
#ifdef DEBUG_OUTPUT
        debug_num_rays++;
#endif
        auto ray = acquireRay();
        ray->setStart(_point_source_locations[src_index], elem);

        ray->setStartingEndPoint(global_spatial_q_points[i]);

        const auto dir = (global_spatial_q_points[i] - _point_source_locations[src_index]).unit();

        for (unsigned int g = 0u; g < _num_groups; ++g)
        {
          ray->data(_source_spatial_weights[g]) =
              global_spatial_q_weights[i] *
              computeSHSource(src_index, _point_source_moments, _point_source_anisotropy, dir, g);
        }

        ray->data(_target_in_element) = -1.0;

        moveRayToBuffer(ray);
      }
    }
  }

  // Surface sources.
  if (_source_boundary_names.size() > 0u)
  {
    const auto & bnd_ids = _mesh.getBoundaryIDs(_source_boundary_names);
    std::unordered_map<unsigned int, std::vector<std::pair<const Elem *, unsigned int>>>
        local_source_bnd_elements;

    for (unsigned int i = 0u; i < _source_boundary_names.size(); ++i)
      local_source_bnd_elements.emplace(i, std::vector<std::pair<const Elem *, unsigned int>>());

    const auto & surface_q_points = _face_fe->get_xyz();
    const auto & surface_q_weights = _face_fe->get_JxW();

    // Find all of the elements we own that contain a surface source.
    for (const auto & b_elem : *_mesh.getBoundaryElementRange())
    {
      for (unsigned int i = 0u; i < bnd_ids.size(); ++i)
      {
        const auto elem = b_elem->_elem;
        const auto side = b_elem->_side;
        if (elem->processor_id() == _pid && b_elem->_bnd_id == bnd_ids[i])
        {
          local_source_bnd_elements[i].emplace_back(std::make_pair(elem, side));
          _face_fe->reinit(elem, side);
          num_surface_source_points += surface_q_points.size();
        }
      }
    }

    total_num_rays += num_surface_source_points * global_spatial_q_points.size();
    reserveRayBuffer(num_surface_source_points * global_spatial_q_points.size());

    // Now that we've found the surface sources we own, we set up the rays.
    for (unsigned int i = 0u; i < global_spatial_q_points.size(); ++i)
    {
      for (const auto & [src_index, elem_vec] : local_source_bnd_elements)
      {
        for (const auto & [elem, side] : elem_vec)
        {
          _face_fe->reinit(elem, side);
          for (unsigned int j = 0u; j < surface_q_points.size(); ++j)
          {
#ifdef DEBUG_OUTPUT
            debug_num_rays++;
#endif
            auto ray = acquireRay();

            ray->setStart(surface_q_points[j], elem);
            ray->setStartingEndPoint(global_spatial_q_points[i]);

            const auto dir = (global_spatial_q_points[i] - surface_q_points[j]).unit();

            for (unsigned int g = 0u; g < _num_groups; ++g)
            {
              ray->data(_source_spatial_weights[g]) =
                  global_spatial_q_weights[i] * surface_q_weights[j] *
                  computeSHSource(
                      src_index, _boundary_source_moments, _boundary_source_anisotropy, dir, g);
            }

            ray->data(_target_in_element) = -1.0;

            moveRayToBuffer(ray);
          }
        }
      }
    }
  }

  // Volume sources.
  if (_volume_source_blocks.size() > 0u)
  {
    const auto elem_ids = _mesh.getSubdomainIDs(_volume_source_blocks);
    std::unordered_map<unsigned int, std::vector<const Elem *>> local_source_elements;

    for (unsigned int i = 0u; i < _volume_source_blocks.size(); ++i)
      local_source_elements.emplace(i, std::vector<const Elem *>());

    const auto & source_q_points = _volume_fe->get_xyz();
    const auto & source_q_weights = _volume_fe->get_JxW();

    // Find all of the elements we own that contain a volumetric source.
    for (const auto & elem : *_mesh.getActiveLocalElementRange())
    {
      for (unsigned int i = 0u; i < elem_ids.size(); ++i)
      {
        if (elem->processor_id() == _pid && elem->subdomain_id() == elem_ids[i])
        {
          local_source_elements[i].emplace_back(elem);
          _volume_fe->reinit(elem);
          num_volume_source_points += source_q_points.size();
        }
      }
    }

    // Subtracting num_volume_src_points * num_volume_src_points to avoid allocating extra rays for
    // the within-element scenario.
    total_num_rays += num_volume_source_points * global_spatial_q_points.size();
    reserveRayBuffer(num_volume_source_points * global_spatial_q_points.size());

    std::size_t num_reallocations = 0u;

    // Now that we've found the volume sources we own, we set up the rays.
    // Out of element contributions go first.
    for (unsigned int i = 0u; i < global_spatial_q_points.size(); ++i)
    {
      for (const auto & [src_index, elem_vec] : local_source_elements)
      {
        for (const auto elem : elem_vec)
        {
          _volume_fe->reinit(elem);
          if (elem->contains_point(global_spatial_q_points[i]))
          {
            num_reallocations += source_q_points.size();
            continue;
          }
          for (unsigned int j = 0u; j < source_q_points.size(); ++j)
          {
#ifdef DEBUG_OUTPUT
            debug_num_rays++;
#endif
            auto ray = acquireRay();

            ray->setStart(source_q_points[j], elem);
            ray->setStartingEndPoint(global_spatial_q_points[i]);

            const auto dir = (global_spatial_q_points[i] - source_q_points[j]).unit();

            for (unsigned int g = 0u; g < _num_groups; ++g)
            {
              ray->data(_source_spatial_weights[g]) =
                  global_spatial_q_weights[i] * source_q_weights[j] *
                  computeSHSource(
                      src_index, _volume_source_moments, _volume_source_anisotropy, dir, g);
            }

            ray->data(_target_in_element) = -1.0;

            moveRayToBuffer(ray);
          }
        }
      }
    }

    // Handle in-element contributions for volumetric sources separately.
    total_num_rays += num_volume_source_points * _num_dir - num_reallocations;
    reserveRayBuffer(num_volume_source_points * _num_dir - num_reallocations);

    for (const auto & [src_index, elem_vec] : local_source_elements)
    {
      for (const auto elem : elem_vec)
      {
        _volume_fe->reinit(elem);
        for (unsigned int i = 0u; i < source_q_points.size(); ++i)
        {
          for (unsigned int n = 0u; n < _num_dir; ++n)
          {
#ifdef DEBUG_OUTPUT
            debug_num_rays++;
#endif
            auto ray = acquireRay();

            RealVectorValue dir;
            if (_dim == 2u)
              dir = RealVectorValue(std::cos(libMesh::pi * _2D_angular_quadrature->root(n)),
                                    std::sin(libMesh::pi * _2D_angular_quadrature->root(n)),
                                    0.0)
                        .unit();
            else
              dir = _3D_angular_quadrature->direction(n).unit();

            const auto angular_weight = _dim == 2u ? libMesh::pi * _2D_angular_quadrature->weight(n)
                                                   : _3D_angular_quadrature->weight(n);

            ray->setStart(source_q_points[i], elem);
            ray->setStartingDirection(dir);

            for (unsigned int g = 0u; g < _num_groups; ++g)
            {
              ray->data(_source_spatial_weights[g]) =
                  angular_weight * source_q_weights[i] *
                  computeSHSource(
                      src_index, _volume_source_moments, _volume_source_anisotropy, dir, g);
            }

            ray->data(_target_in_element) = 1.0;

            moveRayToBuffer(ray);
          }
        }
      }
    }
  }

  _comm.sum(total_num_rays);
  _comm.sum(num_point_source_points);
  _comm.sum(num_surface_source_points);
  _comm.sum(num_volume_source_points);
#ifdef DEBUG_OUTPUT
  _comm.sum(debug_num_rays);
#endif
  _console << "UncollidedFluxRayStudy generated a total of " << total_num_rays << " rays:\n"
           << " - " << global_spatial_q_points.size() << " target points;\n"
           << " - " << num_point_source_points << " point source points;\n"
           << " - " << num_surface_source_points << " surface source points;\n"
           << " - " << num_volume_source_points << " volume source points;\n"
           << " - " << num_volume_source_points * _num_dir << " in-cell volume source rays."
#ifdef DEBUG_OUTPUT
           << "\n - " << debug_num_rays << " rays generated (debug count)"
#endif
           << std::endl;
}
