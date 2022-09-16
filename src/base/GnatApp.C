#include "GnatApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
GnatApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  params.set<bool>("use_legacy_uo_initialization") = false;
  params.set<bool>("use_legacy_uo_aux_computation") = false;
  params.set<bool>("use_legacy_output_syntax") = false;
  params.set<bool>("use_legacy_material_output") = false;

  return params;
}

GnatApp::GnatApp(InputParameters parameters) : MooseApp(parameters)
{
  GnatApp::registerAll(_factory, _action_factory, _syntax);
}

GnatApp::~GnatApp() {}

static void
associateSyntaxInner(Syntax & syntax, ActionFactory & /*action_factory*/)
{
  // Common syntax.
  syntax.registerActionSyntax("CommonGnatAction", "NeutronActivationStudy");

  // Neutron transport syntax.
  syntax.registerActionSyntax("NeutronTransportAction", "NeutronActivationStudy/TransportSystem");

  // Mass transport syntax.
  syntax.registerActionSyntax("SetupIsotopeSystemAction", "NeutronActivationStudy/IsotopeSystem");
  syntax.registerActionSyntax("AddMobileIsotopeAction",
                              "NeutronActivationStudy/IsotopeSystem/AddMobileIsotopes/*");
}

void
GnatApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAll(f, af, syntax);
  Registry::registerObjectsTo(f, {"GnatApp"});
  Registry::registerActionsTo(af, {"GnatApp"});

  /* register custom execute flags, action syntax, etc. here */
  associateSyntaxInner(syntax, af);
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
