//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GnatTestApp.h"
#include "MooseInit.h"
#include "Moose.h"
#include "MooseApp.h"
#include "AppFactory.h"

// Create a performance log
PerfLog Moose::perf_log("Gnat");

// Begin the main program.
int
main(int argc, char * argv[])
{
  // Initialize MPI, solvers and MOOSE
  MooseInit init(argc, argv);

  // Register this application's MooseApp and any it depends on
  GnatTestApp::registerApps();

  // Create an instance of the application and store it in a smart pointer for easy cleanup
  std::shared_ptr<MooseApp> app = AppFactory::createAppShared("GnatTestApp", argc, argv);

  app->_console <<
  "-----------------------------------------------------\n"
  "---------Generic Neutron Activation Toolkit----------\n"
  "-----------------------------------------------------\n"
  "                                               ,----,\n"
  "                     ,--.                    ,/   .`|\n"
  "  ,----..          ,--.'|   ,---,          ,`   .'  :\n"
  " /   /   \\     ,--,:  : |  '  .' \\       ;    ;     /\n"
  "|   :     : ,`--.'`|  ' : /  ;    '.   .'___,/    ,'\n"
  ".   |  ;. / |   :  :  | |:  :       \\  |    :     |\n"
  ".   ; /--`  :   |   \\ | ::  |   /\\   \\ ;    |.';  ;\n"
  ";   | ;  __ |   : '  '; ||  :  ' ;.   :`----'  |  |\n"
  "|   : |.' .''   ' ;.    ;|  |  ;/  \\   \\   '   :  ;\n"
  ".   | '_.' :|   | | \\   |'  :  | \\  \\ ,'   |   |  '\n"
  "'   ; : \\  |'   : |  ; .'|  |  '  '--'     '   :  |\n"
  "'   | '/  .'|   | '`--'  |  :  :           ;   |.'\n"
  "|   :    /  '   : |      |  | ,'           '---'\n"
  " \\   \\ .'   ;   |.'      `--''\n"
  "  `---`     '---'\n"
  "-----------------------------------------------------"
  << std::endl;

  // Execute the application
  app->run();

  return 0;
}
