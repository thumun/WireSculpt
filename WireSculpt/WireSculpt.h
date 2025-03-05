#pragma once

#include <maya/MArgList.h>
#include <maya/MObject.h>
#include <maya/MGlobal.h>
#include <maya/MPxCommand.h>
#include <maya/MSyntax.h>
#include <maya/MArgDatabase.h>
#include <maya/MArgParser.h>


class WireSculpt : public MPxCommand
{
public:
	WireSculpt() {};
	virtual MStatus doIt(const MArgList& args);
	static void* creator();
	static MSyntax newSyntax();
};
