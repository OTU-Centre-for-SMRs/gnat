// Implements Chebyshev polynomials with the intent of computing integrals with
// quadrature sets.
#include "ChebyshevPolynomial.h"

ChebyshevPolynomial::ChebyshevPolynomial(unsigned int degree) : _degree(std::move(degree))
{
  _roots.resize(_degree, 0.0);
  _angular_roots.resize(_degree, 0.0);
  _weights.resize(_degree, libMesh::pi / static_cast<Real>(_degree));

  for (unsigned int i = 1; i <= _degree; ++i)
  {
    _angular_roots[i - 1u] =
        (2 * static_cast<Real>(i) - 1.0) / (2.0 * static_cast<Real>(_degree)) * libMesh::pi;
    _roots[i - 1u] = std::cos(_angular_roots[i - 1u]);
  }
}
