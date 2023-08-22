#pragma once

// Major axis for quadrature sets.
enum class MajorAxis
{
  X = 0u,
  Y = 1u,
  Z = 2u
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

// An enum for particle types. Either neutron or photon.
enum class Particletype
{
  Neutron = 0u,
  GammaPhoton = 1u
};
// enum class Particletype

// An enum for radiation transport schemes.
enum class TransportScheme
{
  SAAFCFEM = 0u,
  DiffusionApprox = 1u,
  FluxMomentTransfer = 2u
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

// An enum to make the parsing code for cross-section files cleaner.
enum class PropertyType
{
  InvVelocity = 0u,
  SigmaR = 1u,
  SigmaA = 2u,
  SigmaS = 3u,
  SigmaSMatrix = 4u,
  NuSigmaF = 5u,
  ChiF = 6u
}; // enum class PropertyType

// An enum to make the parsing code for cross-section files cleaner.
enum class CrossSectionSource
{
  Detect = 0u,
  Gnat = 1u,
  OpenMC = 2u
}; // enum class CrossSectionSource

// An enum to indicate if cross-sections are macroscopic or microscopic.
enum class CrossSectionType
{
  Micro = 0u,
  Macro = 1u
}; // enum class CrossSectionType
