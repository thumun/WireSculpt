#include "WireSculpt.h"
#include "WireSculptPlugin.h"
#include <maya/MFnPlugin.h>
#include "WireSculptNode.h"

// define EXPORT for exporting dll functions
#define EXPORT _declspec(dllexport)

const char* fileFlag = "-fp";
const char* fileLongFlag = "-filepath";

const char* attractStrengthFlag = "-asr";
const char* attractStrengthLongFlag = "-attractStrength";

const char* attractSteepFlag = "-ase";
const char* attractSteepLongFlag = "-attractSteep";

const char* repulsionStrengthFlag = "-rsr";
const char* repulsionStrengthLongFlag = "-repulsionStrength";

const char* repulsionSteepFlag = "-rse";
const char* repulsionSteepLongFlag = "-repulsionSteep";

const char* expandFlag = "-e";
const char* expandLongFlag = "-expand";

const char* keepFlag = "-k";
const char* keepLongFlag = "-keep";

const char* lambdaFlag = "-l";
const char* lambdaLongFlag = "-lambda";

const char* thicknessFlag = "-t";
const char* thicknessLongFlag = "-thickness";

void* WireSculpt::creator()
{
	return new WireSculpt;
}

MSyntax WireSculpt::newSyntax()
{
	MSyntax syntax;

	syntax.addFlag(fileFlag, fileLongFlag, MSyntax::kString);
	syntax.addFlag(attractStrengthFlag, attractStrengthLongFlag, MSyntax::kDouble);
	syntax.addFlag(attractSteepFlag, attractSteepLongFlag, MSyntax::kDouble);
	syntax.addFlag(repulsionStrengthFlag, repulsionStrengthLongFlag, MSyntax::kDouble);
	syntax.addFlag(repulsionSteepFlag, repulsionSteepLongFlag, MSyntax::kDouble);
	syntax.addFlag(expandFlag, expandLongFlag, MSyntax::kLong);
	syntax.addFlag(keepFlag, keepLongFlag, MSyntax::kLong);
	syntax.addFlag(lambdaFlag, lambdaLongFlag, MSyntax::kDouble);
	syntax.addFlag(thicknessFlag, thicknessLongFlag, MSyntax::kDouble);

	return syntax;
}

// Plugin doIt function
MStatus WireSculpt::doIt(const MArgList& argList)
{
	MStatus status;

	MString file = "";
	double attractStrength = 1.0;
	double attractSteep = 0.0;
	double repulStrength = 1.0;
	double repulSteep = 0.0;
	int expand = 0;
	int keep = 0;
	double lambda = 1.0;
	double thickness = 1.0;


	MArgParser argData(syntax(), argList, &status);

	if (argData.isFlagSet(fileFlag)) {
		argData.getFlagArgument(fileFlag, 0, file);
	}
	else {
		MGlobal::displayInfo("Flag -fp is NOT set.");
	}

	if (argData.isFlagSet(attractStrengthFlag)) {
		argData.getFlagArgument(attractStrengthFlag, 0, attractStrength);
	}
	else {
		MGlobal::displayInfo("Flag -asr is NOT set.");
	}

	if (argData.isFlagSet(attractSteepFlag)) {
		argData.getFlagArgument(attractSteepFlag, 0, attractSteep);
	}
	else {
		MGlobal::displayInfo("Flag -ase is NOT set.");
	}

	if (argData.isFlagSet(repulsionStrengthFlag)) {
		argData.getFlagArgument(repulsionStrengthFlag, 0, repulStrength);
	}
	else {
		MGlobal::displayInfo("Flag -rsr is NOT set.");
	}

	if (argData.isFlagSet(repulsionSteepFlag)) {
		argData.getFlagArgument(repulsionSteepFlag, 0, repulSteep);
	}
	else {
		MGlobal::displayInfo("Flag -rse is NOT set.");
	}

	if (argData.isFlagSet(expandFlag)) {
		argData.getFlagArgument(expandFlag, 0, expand);
	}
	else {
		MGlobal::displayInfo("Flag -e is NOT set.");
	}

	if (argData.isFlagSet(keepFlag)) {
		argData.getFlagArgument(keepFlag, 0, keep);
	}
	else {
		MGlobal::displayInfo("Flag -k is NOT set.");
	}

	if (argData.isFlagSet(lambdaFlag)) {
		argData.getFlagArgument(lambdaFlag, 0, lambda);
	}
	else {
		MGlobal::displayInfo("Flag -l is NOT set.");
	}

	if (argData.isFlagSet(thicknessFlag)) {
		argData.getFlagArgument(thicknessFlag, 0, thickness);
	}
	else {
		MGlobal::displayInfo("Flag -t is NOT set.");
	}

	// flag values 
	MString print = "filepath: ";
	print += file;
	print += "\n attract stren: ";
	print += attractStrength;
	print += "\n attract steep: ";
	print += attractSteep;
	print += "\n repulsion stren: ";
	print += repulStrength;
	print += "\n repulsion steep: ";
	print += repulSteep;
	print += "\n expand: ";
	print += expand;
	print += "\n keep: ";
	print += keep;
	print += "\n lambda: ";
	print += lambda;
	print += "\n thickness: ";
	print += thickness;
	MGlobal::displayInfo(print);

	// process file
	WireSculptPlugin ws = WireSculptPlugin();
	bool returnVal = ws.ProcessFile(file.asChar());
	MString errorMsg = "\n";
	if (returnVal) {
		errorMsg += "file processed successfully";
		MGlobal::displayInfo(errorMsg);
		MGlobal::displayInfo(MString() + ws.GetVerticies()->size());
	}
	else {
		errorMsg += "issue with file format or contents";
		MGlobal::displayInfo(errorMsg);
	}

	return MStatus::kSuccess;
}

// creates menu to open WireSculpt UI
void createMenu(MString& melPath) {
	MString menuName = "WireSculptMenu";

	MGlobal::executeCommand(MString("if (`menu -exists ") + menuName + MString("`) deleteUI ") + menuName + MString(";"));

	MGlobal::executeCommand(MString("menu -parent MayaWindow -label \"WireSculpt\" ") + menuName + MString(";"));

	MString sourceCmd = "source \"" + melPath + "\"";
	MStatus sourceStatus = MGlobal::executeCommand(sourceCmd);

	if (sourceStatus != MS::kSuccess) {
		MGlobal::displayError("[WireSculpt] Failed to source MEL script: " + melPath);
		return;
	}

	MGlobal::executeCommand(
		"menuItem -parent " + menuName +
		" -label \"Open WireSculpt\" " +
		"-command \"CreateWireSculpt;\""
	);
}

// Initialize Maya Plugin upon loading
EXPORT MStatus initializePlugin(MObject obj)
{
	MStatus status;
	MFnPlugin plugin(obj, "CIS660", "1.0", "Any");

	status = plugin.registerCommand("WireSculpt", WireSculpt::creator, WireSculpt::newSyntax);
	if (!status)
		status.perror("registerCommand failed");

	// Run with the MEL command createNode WireSculptNode
	status = plugin.registerNode("WireSculptNode", WireSculptNode::id,
		WireSculptNode::creator, WireSculptNode::initialize);
	if (!status) {
		status.perror("registerNode");
		return status;
	}

	MString pluginPath = plugin.loadPath();
	MString melPath = pluginPath + "/wireUI.mel";

	createMenu(melPath);

	return status;
}
// Cleanup Plugin upon unloading
EXPORT MStatus uninitializePlugin(MObject obj)
{
	MStatus status;
	MFnPlugin plugin(obj);

	status = plugin.deregisterNode(WireSculptNode::id);
	if (!status) {
		status.perror("deregisterNode");
		return status;
	}

	MGlobal::executeCommand("if (`menu -exists WireSculptMenu`) deleteUI WireSculptMenu;");
	MGlobal::executeCommand("if (`menu -exists wireSculptMenuUI`) deleteUI wireSculptMenuUI;");

	status = plugin.deregisterCommand("WireSculpt");
	if (!status)
		status.perror("deregisterCommand failed");
	return status;
}