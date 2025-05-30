//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GnatApp.h"
#include "MooseMain.h"

// Create a performance log
PerfLog Moose::perf_log("Gnat");

// Begin the main program.
int
main(int argc, char * argv[])
{
  Moose::main<GnatApp>(argc, argv);

  return 0;
}
