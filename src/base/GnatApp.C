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
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;

  return params;
}

GnatApp::GnatApp(InputParameters parameters) : MooseApp(parameters)
{
  GnatApp::registerAll(_factory, _action_factory, _syntax);

  _console << "-----------------------------------------------------\n"
              "------General Neutral Particle Analysis Toolkit------\n"
              "-----------------------------------------------------\n"
              "                     ,--.\n"
              "  ,----..          ,--.'|                      ,----,\n"
              " /   /   \\     ,--,:  : |   ,---,            ,/   .`|\n"
              "|   :     : ,`--.'`|  ' :  '  .' \\         ,`   .'  :\n"
              ".   |  ;. / |   :  :  | | /  ;    '.     ;    ;     /\n"
              ".   ; /--`  :   |   \\ | ::  :   .   \\  .'___,/    ,'\n"
              ";   | ;  __ |   : '  '; |:  |  / \\   \\ |    :     |\n"
              "|   : |.' .''   ' ;.    ;|  :  |--\\   :;    |.';  ;\n"
              ".   | '_.' :|   | | \\   ||  |  /\\  \\   `----'  |  |\n"
              "'   ; : \\  |'   : |  ; .''  :  | \\  \\ ,'   '   :  ;\n"
              "'   | '/  .'|   | '`--'  |  |  '  '--'     |   |  '\n"
              "|   :    /  '   : |      |  :  :           '   :  |\n"
              " \\   \\ .'   ;   |.'      |  | ,'           ;   |.'\n"
              "  `---`     '---'        `--''             '---'\n"
              "-----------------------------------------------------"
           << std::endl;
}

GnatApp::~GnatApp() {}

static void
associateSyntaxInner(Syntax & syntax, ActionFactory & /*action_factory*/)
{
  // TransportSystem syntax.
  syntax.registerActionSyntax("AddTransportMaterialAction", "TransportMaterials/*");
  syntax.registerActionSyntax("TransportAction", "TransportSystems/*");
  syntax.registerActionSyntax("UncollidedFluxAction", "UncollidedFlux/*");

  syntax.registerActionSyntax("CardinalTransportMaterialAction", "CardinalMGXS");

  // Depletion library.
  syntax.registerActionSyntax("DepletionLibraryAction", "DepletionLibrary");

  // Mobile depletion syntax.
  syntax.registerActionSyntax("MobileDepletionSystemAction", "MobileDepletionSystem");
  syntax.registerActionSyntax("TracerDepletionSystemAction", "TracerDepletionSystem");
}

void
GnatApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAllObjects<GnatApp>(f, af, syntax);
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
