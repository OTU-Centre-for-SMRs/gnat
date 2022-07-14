#pragma once

// Major axis for quadrature sets.
enum class MajorAxis
{
  X = 0,
  Y = 1,
  Z = 2
}; // enum class MajorAxis

// Problem type to determine the appropriate quadrature set and spherical
// harmonics expansion.
// TODO: Cylindrical and spherical coordiante systems + meshes.
enum class ProblemType
{
  Cartesian1D = 0u,
  Cartesian2D = 1u,
  Cartesian3D = 2u
}; // enum class ProblemType
