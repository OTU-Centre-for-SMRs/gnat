#pragma once

#include "ADKernel.h"

// This kernel evaluates the full scattering contribution. This way the complete
// matrix is assembled and source iteration is avoided. This kernel should only
// be used for debugging purposes as a fully assembled matrix solve for the transport
// equation is quite slow.
// TODO: Finish this kernel.
class ADNeutronScattering : public ADKernel
{
public:

protected:

};
