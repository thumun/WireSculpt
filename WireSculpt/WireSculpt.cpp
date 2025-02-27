#include "WireSculpt.h"
#include <maya/MFnPlugin.h>
#include "WireSculptNode.h"

// define EXPORT for exporting dll functions
#define EXPORT _declspec(dllexport)

const char* nameFlag = "-n";
const char* nameLongFlag = "-name";
const char* idFlag = "-id";
const char* idLongFlag = "-idnumber";

void* WireSculpt::creator()
{
	return new WireSculpt;
}

MSyntax WireSculpt::newSyntax()
{
	MSyntax syntax;
	syntax.addFlag(nameFlag, nameLongFlag, MSyntax::kString);
	syntax.addFlag(idFlag, idLongFlag, MSyntax::kUnsigned);
	return syntax;
}

// Plugin doIt function
MStatus WireSculpt::doIt(const MArgList& argList)
{
	MGlobal::displayInfo("about to add script");

	MStatus status;

	MString firstFlagArg;
	unsigned int secondFlagArg{};

	MArgParser argData(syntax(), argList, &status);

	if (argData.isFlagSet(nameFlag)) {
		argData.getFlagArgument(nameFlag, 0, firstFlagArg);
	}
	else {
		MGlobal::displayInfo("Flag -name is NOT set.");
	}

	if (argData.isFlagSet(idFlag)) {
		argData.getFlagArgument(idFlag, 0, secondFlagArg);
	}
	else {
		MGlobal::displayInfo("Flag -id is NOT set.");
	}

	MGlobal::displayInfo("pop-up menu");

	MString result;

	MString dialog = "confirmDialog -title \"Hello Maya\" -message \"Name: ";
	dialog += firstFlagArg;
	dialog += "\\nID: ";
	dialog += secondFlagArg;
	dialog += "\" -button \"OK\" -defaultButton \"OK\" -dismissString \"Cancel\";";
	MGlobal::executeCommand(dialog, result);

	return MStatus::kSuccess;
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

	status = plugin.deregisterCommand("WireSculpt");
	if (!status)
		status.perror("deregisterCommand failed");
	return status;
}