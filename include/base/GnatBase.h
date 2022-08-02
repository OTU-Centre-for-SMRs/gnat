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

enum class Scheme
{
  SAAFCFEM = 0u,
  UpwindingDFEM = 1u
}; // enum class Scheme

enum class ExecutionType
{
  SteadySource = 0u,
  Transient = 1u
}; // enum class ExecutionType

enum class DebugVerbosity
{
  Level0 = 0u,
  Level1 = 1u
}; // Enum class DebugVerbosity
