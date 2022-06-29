#include "visitorextract.h"
#include "iso_method_ours.h"
#include "tet_arrays.h"

TraversalData::TraversalData(TNode *t)
	{
		n = t;
		depth = 0;
	}

	
void TraversalData::gen_trav(TraversalData &c, Index i)
	{
		if (!n->is_leaf())
		{
			c.n = n->children[i];
			c.depth = depth+1;
		}
		else
		{
			c = *this;
		}
	}

	VisitorExtract::VisitorExtract(Mesh *m_)
	{
		m = m_;
	}

	bool VisitorExtract::on_node(TraversalData &td){return !td.n->is_leaf();}

	bool VisitorExtract::on_edge(TraversalData &td00, TraversalData &td10, TraversalData &td01, TraversalData &td11, char orient)
	{
		if (td00.n->is_leaf() && td10.n->is_leaf() && td01.n->is_leaf() && td11.n->is_leaf())
		{
			// determine which node is smallest (ordered bitwise)
			TNode *n[4] = {td00.n, td10.n, td01.n, td11.n};
			int depths[4] = {td00.depth, td10.depth, td01.depth, td11.depth};
			int small = 0;
			for (int i = 1; i < 4; i++)
			{
				if (depths[i] > depths[small])
					small = i;
			}

			// calc edge vertex
			int edge_idx = cube_orient2edge[orient][small ^ 3];

			vect4d &edgeMid = n[small]->edges[edge_idx];
			vect4d &edgeLow = n[small]->verts[cube_edge2vert[edge_idx][0]];
			vect4d &edgeHigh = n[small]->verts[cube_edge2vert[edge_idx][1]];

			// init common tetrahedron points
			vect4d p[6];
			p[0] = edgeLow; 
			p[1] = edgeMid;
			p[4] = edgeMid; 
			p[5] = edgeHigh;

#ifdef JOIN_VERTS
			vect3d topoLow = n[small]->verts[cube_edge2vert[edge_idx][0]];
			vect3d topoHigh = n[small]->verts[cube_edge2vert[edge_idx][1]];
			vect3d topoMid = (topoLow + topoHigh) * .5;

			vect3d topo[6];
			topo[0] = topoLow; 
			topo[1] = topoMid;
			topo[4] = topoMid; 
			topo[5] = topoHigh;
#endif

			// create pyramids in each cell
			for (int i = 0; i < 4; i++)
			{
				if (n[i]->is_outside())
					continue;

				static int faceTable[3][4][2] =
				{
					{
						{1,0},{4,0},{1,5},{4,5}
					},
					{
						{2,0},{3,0},{2,5},{3,5}
					},
					{
						{2,1},{3,1},{2,4},{3,4}
					}
				};

				static bool flipTable[3][4] =
				{
					{true,false,false,true},
					{false,true,true,false},
					{true,false,false,true}
				};

				// find min size face
				const int i1 = i ^ 1;
				const int i2 = i ^ 2;
				bool do1 = n[i] != n[i1];
				bool do2 = n[i] != n[i2];
				vect4d face1, face2;

				if (depths[i] == depths[i1])
					face1 = (n[i]->faces[faceTable[orient][i^3][0]] + n[i1]->faces[cube_face2opposite[faceTable[orient][i^3][0]]]) * .5;
				else if (depths[i] > depths[i1])
					face1 = n[i]->faces[faceTable[orient][i^3][0]];
				else
					face1 = n[i1]->faces[cube_face2opposite[faceTable[orient][i^3][0]]];

				if (depths[i] == depths[i2])
					face2 = (n[i]->faces[faceTable[orient][i^3][1]] + n[i2]->faces[cube_face2opposite[faceTable[orient][i^3][1]]]) * .5;
				else if (depths[i] > depths[i2])
					face2 = n[i]->faces[faceTable[orient][i^3][1]];
				else
					face2 = n[i2]->faces[cube_face2opposite[faceTable[orient][i^3][1]]];

#ifdef JOIN_VERTS
				vect3d topoFace1, topoFace2;

				if (depths[i] > depths[i1])
				{
					int face = faceTable[orient][i^3][0];
					topoFace1 = (n[i]->verts[cube_face2vert[face][0]] + n[i]->verts[cube_face2vert[face][2]]) * .5;
				}
				else
				{
					int face = cube_face2opposite[faceTable[orient][i^3][0]];
					topoFace1 = (n[i1]->verts[cube_face2vert[face][0]] + n[i1]->verts[cube_face2vert[face][2]]) * .5;
				}

				if (depths[i] > depths[i2])
				{
					int face = faceTable[orient][i^3][1];
					topoFace2 = (n[i]->verts[cube_face2vert[face][0]] + n[i]->verts[cube_face2vert[face][2]]) * .5;
				}
				else
				{
					int face = cube_face2opposite[faceTable[orient][i^3][1]];
					topoFace2 = (n[i2]->verts[cube_face2vert[face][0]] + n[i2]->verts[cube_face2vert[face][2]]) * .5;
				}
#endif

				if (flipTable[orient][i])
				{
					if (do1)
					{
						p[2] = face1; p[3] = n[i]->node;		
#ifdef JOIN_VERTS
						topo[2] = topoFace1; topo[3] = n[i]->node;
						processTet(p, topo);
						processTet(p + 2, topo + 2);
#else
						processTet(p);
						processTet(p + 2);
#endif
					}
					if (do2)
					{
						p[2] = n[i]->node; p[3] = face2;	
#ifdef JOIN_VERTS
						topo[2] = n[i]->node; topo[3] = topoFace2;
						processTet(p, topo);
						processTet(p + 2, topo + 2);
#else
						processTet(p);
						processTet(p + 2);
#endif
					}
				}
				else
				{
					if (do1)
					{
						p[2] = n[i]->node; p[3] = face1;	
#ifdef JOIN_VERTS
						topo[2] = n[i]->node; topo[3] = topoFace1;
						processTet(p, topo);
						processTet(p + 2, topo + 2);
#else
						processTet(p);
						processTet(p + 2);
#endif
					}
					if (do2)
					{
						p[2] = face2; p[3] = n[i]->node;
#ifdef JOIN_VERTS
						topo[2] = topoFace2; topo[3] = n[i]->node;
						processTet(p, topo);
						processTet(p + 2, topo + 2);
#else
						processTet(p);
						processTet(p + 2);
#endif
					}
				}

			}
			return false;
		}

		return true;
	}

	bool VisitorExtract::on_face(TraversalData &td0, TraversalData &td1, char orient)
	{
		return !(td0.n->is_leaf() && td1.n->is_leaf());
	}
#ifdef JOIN_VERTS
	void VisitorExtract::processTet(vect4d *p, vect3d *topo)
#else
	void VisitorExtract::processTet(vect4d *p)
#endif
	{
		// determine index in lookup table
		int idx = (p[0].v[3] >= 0 ? 1 : 0) | 
				  (p[1].v[3] >= 0 ? 2 : 0) | 
				  (p[2].v[3] >= 0 ? 4 : 0) | 
				  (p[3].v[3] >= 0 ? 8 : 0);

		for (int i = 0; tet_tris[idx][i] != -1; i += 3)
		{
			vect< 3, vect3d > verts;
#ifdef JOIN_VERTS
			vect< 3, TopoEdge > topoEdges;
#endif
			for (int j = 0; j < 3; j++)
			{
				int edge = tet_tris[idx][i+j];
				verts[j] = findZero(
					p[ tet_edge2vert[edge][0] ], 
					p[ tet_edge2vert[edge][1] ]);
				
#ifdef JOIN_VERTS
				topoEdges[j].v[0] = topo[ tet_edge2vert[edge][0] ];
				topoEdges[j].v[1] = topo[ tet_edge2vert[edge][1] ];
				topoEdges[j].fix();
#endif
			}

			m->tris.push_back(verts);
#ifdef JOIN_VERTS
			m->topoTris.push_back(topoEdges);
#endif
		}
	}
