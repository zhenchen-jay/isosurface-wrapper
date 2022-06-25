#pragma once 
#include "iso_common.h"
#include "vect.h"
#include "index.h"

struct Mesh;
struct TNode;

struct TraversalData
{
	TraversalData(){}
	TraversalData(TNode *t);

	TNode *n; // node
	char depth; 
	
	void gen_trav(TraversalData &c, Index i);
};

#ifdef JOIN_VERTS
#include <algorithm>

using namespace std;

struct TopoEdge
{
	TopoEdge(){}
	TopoEdge(vect3f &c, vect3f &d)
	{
		v[0] = c;
		v[1] = d;
		fix();
	}

	vect3f v[2];

	void fix()
	{
		if (v[1] < v[0])
		{
			swap(v[0], v[1]);
		}
	}
	bool operator<(const TopoEdge &a) const 
	{
		return lexicographical_compare(v, v+2, a.v, a.v+2);
	}
};
#endif

struct VisitorExtract
{
	VisitorExtract(Mesh *m_);
	Mesh *m;

	bool on_vert(TraversalData &a, TraversalData &b, TraversalData &c, TraversalData &d,
		TraversalData &aa, TraversalData &ba, TraversalData &ca, TraversalData &da){return false;}

	bool on_node(TraversalData &td);

	bool on_edge(TraversalData &td00, TraversalData &td10, TraversalData &td01, TraversalData &td11, char orient);

	bool on_face(TraversalData &td0, TraversalData &td1, char orient);

#ifdef JOIN_VERTS
	void processTet(vect4f *p, vect3f *topo);
#else
	void processTet(vect4f *p);
#endif
};