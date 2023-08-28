#pragma once

#include "IntegralRayKernelBase.h"

#include "MooseVariableFE.h"
#include "MooseVariableInterface.h"

class AuxiliarySystem;

/*
 * A class which computes the uncollided flux for point/volume/surface sources using raytracing.
 * This has to combine functionality from integral ray kernels (computing the optical depth) and
 * residual ray kernels (computing the uncollided flux).

 * This should really be an auxkernel, but MOOSE currently doesn't support higher order and
 * non-nodal ray auxkernels.
 */
class UncollidedFluxRayKernel : public IntegralRayKernelBase,
                                public MooseVariableInterface<RealEigenVector>
{
public:
  static InputParameters validParams();

  UncollidedFluxRayKernel(const InputParameters & parameters);

  /**
   * Gets the name of the Ray data associated with the integral accumulated by this RayKernel
   */
  std::string integralRayDataName() const { return _name + "_value"; }

  /**
   * Gets the variable this kernel operates on.
   */
  MooseVariableFE<RealEigenVector> & variable() { return _var; }

  void addValue(const RealEigenVector & value);

  void onSegment() override final;

protected:
  static void cartesianToSpherical(const RealVectorValue & direction, Real & mu, Real & omega);

  // Function to compute the segment contribution to the optical depth.
  void computeSegmentOpticalDepth();

  // Functions to compute the uncollided flux when the ray hits the target point. Two cases:
  // - The source element is equal to the target element, and so we need to remove the Green's
  // function to avoid explosions up to infinity.
  // - The source is not in the target element and so we can use the Green's function.
  void computeUncollidedFluxSourceIsTarget();
  void computeUncollidedFluxSourceNotTarget();

  // The aux system
  AuxiliarySystem & _aux;

  // The AuxVariable this AuxRayKernel contributes to
  MooseVariableFE<RealEigenVector> & _var;

  // The index into the data on the ray that this integral accumulates into.
  std::vector<RayDataIndex> _integral_data_indices;
  // The index into the data on the ray which contains the source intensity and spatial quadratures,
  // pre-multiplied.
  std::vector<RayDataIndex> _source_spatial_weights;
  // The index into the ray which contains a sign to indicate if the target point and source point
  // are in the same element.
  const RayDataIndex _target_in_element;

  // Data required for the optical depth.
  const unsigned int _num_groups;
  const unsigned int _max_eval_anisotropy;
  const unsigned int _num_group_moments;
  // The total cross-section.
  const ADMaterialProperty<std::vector<Real>> & _sigma_t_g;

private:
  /// Spin mutex object for adding values
  static Threads::spin_mutex _add_value_mutex;
}; // class UncollidedFluxRayKernel
