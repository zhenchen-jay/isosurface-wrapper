#include "global.h"
#include "parse.h"
#include "iso_method_ours.h"

Global g;

float pi = 3.1415926535f;

Global::Global()
{
	method = 0;
	ourRoot = 0;
}

