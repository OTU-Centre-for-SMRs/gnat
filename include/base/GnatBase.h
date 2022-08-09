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
// TODO: Cylindrical and spherical coordiante systems.
enum class ProblemType
{
  Cartesian1D = 0u,
  Cartesian2D = 1u,
  Cartesian3D = 2u
}; // enum class ProblemType

// An enum for schemes. Either SAAF with CGFEM or Upwinding with DGFEM.
enum class Scheme
{
  SAAFCFEM = 0u,
  UpwindingDFEM = 1u
}; // enum class Scheme

// An enum for the execution type of the problem. Either steady-state or
// transient.
enum class ExecutionType
{
  SteadySource = 0u,
  Transient = 1u
}; // enum class ExecutionType

// An enum to determine the level of debug output verbosity. Level0 is fully
// verbose. Level1 is minimal.
enum class DebugVerbosity
{
  Level0 = 0u,
  Level1 = 1u
}; // enum class DebugVerbosity
