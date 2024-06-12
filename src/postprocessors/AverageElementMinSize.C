#include "AverageElementMinSize.h"

registerMooseObject("MooseApp", AverageElementMinSize);

InputParameters
AverageElementMinSize::validParams()
{
  InputParameters params = ElementPostprocessor::validParams();
  params.addClassDescription("Computes the average minimum element size.");
  return params;
}

AverageElementMinSize::AverageElementMinSize(const InputParameters & parameters)
  : ElementPostprocessor(parameters)
{
}

void
AverageElementMinSize::initialize()
{
  _total_size = 0;
  _elems = 0;
}

void
AverageElementMinSize::execute()
{
  _total_size += _current_elem->hmin();
  _elems++;
}

Real
AverageElementMinSize::getValue() const
{
  return _total_size / _elems;
}

void
AverageElementMinSize::threadJoin(const UserObject & y)
{
  const auto & pps = static_cast<const AverageElementMinSize &>(y);
  _total_size += pps._total_size;
  _elems += pps._elems;
}

void
AverageElementMinSize::finalize()
{
  gatherSum(_total_size);
  gatherSum(_elems);
}
