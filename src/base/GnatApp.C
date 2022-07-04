#include "GnatApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
GnatApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  return params;
}

GnatApp::GnatApp(InputParameters parameters) : MooseApp(parameters)
{
  GnatApp::registerAll(_factory, _action_factory, _syntax);
}

GnatApp::~GnatApp() { }

void
GnatApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAll(f, af, syntax);
  Registry::registerObjectsTo(f, {"GnatApp"});
  Registry::registerActionsTo(af, {"GnatApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
GnatApp::registerApps()
{
  registerApp(GnatApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
GnatApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  GnatApp::registerAll(f, af, s);
}
extern "C" void
GnatApp__registerApps()
{
  GnatApp::registerApps();
}
