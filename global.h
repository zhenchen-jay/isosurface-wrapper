#pragma once

#include <string>
#include <vector>
#include <math.h>
#include "iso_common.h"
#include "player.h"
#include "vect.h"
#include "visitorextract.h"
using namespace std;


struct SDL_Surface;
struct TerrainWidget;
struct TNode;
struct Widget;


struct Mesh
{
	vector< vect<3, vect3f> > tris;
#ifdef JOIN_VERTS
	vector< vect<3, TopoEdge> > topoTris;
#endif
	vector<vect3f> norms;
};

struct Global
{
	Global();

	TNode *ourRoot;

	Player player, camera;


	int method;
	Mesh ourMesh, rootsMesh, dmcMesh, dcMesh;
};

extern Global g;
