#include "iso_method_ours.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "index.h"
#include "matrix.h"
#include <set>
#include "timer.h"
#include "csgtree.h"
#include "visitorextract.h"
#include "rootfind.h"
#include "traverse.h"

int tree_cells = 0;

void gen_iso_ours()
{
	double t_start = get_time();

	TNode *guide = g.ourRoot;
	TNode *root;
	Mesh *m = &g.ourMesh;

	if (guide)
	{
		printf("-= Our Method =-\n");

		root = guide;
	}
	else
	{
		printf("-= Calculating Tree Structure =-\n");
		root = new TNode();
		g.ourRoot = root;
	}

	for (Index i; i < 8; i++)
		root->verts[i.v](int(i.x)*4-1, int(i.y)*4-1, int(i.z)*4-1);

	vect3d grad[8];
	for (Index i; i < 8; i++)
		csg_root->eval(*(vect3d*)&root->verts[i], root->verts[i].v[3], grad[i]);
	
	root->eval(grad, guide);

	double t_finish = get_time();
	printf("Time generating tree = %f\n", t_finish - t_start);

	// extract surface
	m->tris.reserve(1000000);
	VisitorExtract v(m);
	TraversalData td(root);
	traverse_node<trav_edge>(v, td);

	double t_alldone = get_time();
	printf("Time generating polygons = %f\n", t_alldone - t_finish);
	printf("Time total = %f\n", t_alldone - t_start);

	// calc normals
	m->norms.resize(m->tris.size());
	for (int i = 0; i < m->tris.size(); i++)
	{
		m->norms[i] = (m->tris[i][1] - m->tris[i][0]) % (m->tris[i][2] - m->tris[i][0]);
		m->norms[i].normalize();
	}


	//if (guide)
	//{
	//	// extract surface
	//	m->tris.reserve(1000000);
	//	VisitorExtract v(m);
	//	TraversalData td(root);
	//	traverse_node<trav_edge>(v, td);
	//	
	//	double t_alldone = get_time();
	//	printf("Time generating polygons = %f\n", t_alldone - t_finish);
	//	printf("Time total = %f\n", t_alldone - t_start);

	//	// calc normals
	//	m->norms.resize(m->tris.size());
	//	for (int i = 0; i < m->tris.size(); i++)
	//	{
	//		m->norms[i] = (m->tris[i][1] - m->tris[i][0]) % (m->tris[i][2] - m->tris[i][0]);
	//		m->norms[i].normalize();
	//	}
	//}
	//else
	//{
	//	printf("Cells in tree = %d\n", tree_cells);
	//}
}

bool TNode::changesSign(vect4d *verts, vect4d *edges, vect4d *faces, vect4d &node)
{
	return  sign(verts[0]) != sign(verts[1]) || 
			sign(verts[0]) != sign(verts[2]) || 
			sign(verts[0]) != sign(verts[3]) ||
			sign(verts[0]) != sign(verts[4]) || 
			sign(verts[0]) != sign(verts[5]) || 
			sign(verts[0]) != sign(verts[6]) || 
			sign(verts[0]) != sign(verts[7]) || 

			sign(verts[0]) != sign(faces[1]) || 
			sign(verts[0]) != sign(faces[2]) || 
			sign(verts[0]) != sign(faces[3]) ||
			sign(verts[0]) != sign(faces[4]) || 
			sign(verts[0]) != sign(faces[5]) || 

			sign(verts[0]) != sign(edges[1]) || 
			sign(verts[0]) != sign(edges[2]) || 
			sign(verts[0]) != sign(edges[3]) ||
			sign(verts[0]) != sign(edges[4]) || 
			sign(verts[0]) != sign(edges[5]) || 
			sign(verts[0]) != sign(edges[6]) || 
			sign(verts[0]) != sign(edges[7]) || 
			sign(verts[0]) != sign(edges[8]) || 
			sign(verts[0]) != sign(edges[9]) ||
			sign(verts[0]) != sign(edges[10]) || 
			sign(verts[0]) != sign(edges[11]) || 

			sign(verts[0]) != sign(node);
}

void TNode::eval(vect3d *grad, TNode *guide)
{
	double qef_error = 0;
	if (std::abs(verts[0][0] - 0.5) < 1e-6 && std::abs(verts[0][1] - 0.25) < 1e-6 && std::abs(verts[0][2] - 0.5) < 1e-6
		&& std::abs(verts[7][0] - 0.75) < 1e-6 && std::abs(verts[7][1] - 0.5) < 1e-6 && std::abs(verts[7][2] - 0.75) < 1e-6)
		std::cout << "debug" << std::endl;

	if (!guide || (guide && guide->children[0] == 0))
	{
		tree_cells++;

		// evaluate QEF samples
		vertNode(node, grad, qef_error);
		for (int i = 0; i < 12; i++)
			vertEdge(edges[i], grad, qef_error, i);
		for (int i = 0; i < 6; i++)
			vertFace(faces[i], grad, qef_error, i);
	}

	bool recur;

	if (!guide)
	{
		double cellsize = verts[7][0] - verts[0][0];

		if (is_outside((vect3d&)verts[0], (vect3d&)verts[7]))
			return;

		// check max/min sizes of cells
		static double minsize = pow(.5, DEPTH_MAX);
		bool issmall = cellsize <= minsize;
		if (issmall)
		{
			// it's a leaf
			//printf("leaf size = %f\n", cellsize);
			return;
		}

		static double maxsize = pow(.5, DEPTH_MIN);
		bool isbig = cellsize > maxsize;

		// check for a sign change
		bool signchange = false;
		if (!isbig)
			signchange = changesSign(verts, edges, faces, node);

		// check for qef error
		bool badqef = qef_error/cellsize > 1e-3;

		recur = isbig || (signchange && badqef);
	}
	else
	{
		recur = guide->children[0] != 0;
	}






	if (recur)
	{
		// find points and function values in the subdivided cell
		vect4d p[3][3][3];
		vect3d g[3][3][3];
		for (int x = 0; x < 3; x++)
		{
			for (int y = 0; y < 3; y++)
			{
				for (int z = 0; z < 3; z++)
				{
					p[x][y][z] =(verts[Index(0,0,0)]*(2-x)*(2-y)*(2-z) +
								 verts[Index(0,0,1)]*(2-x)*(2-y)*(z) +
								 verts[Index(0,1,0)]*(2-x)*(y)*(2-z) +
								 verts[Index(0,1,1)]*(2-x)*(y)*(z) +
								 verts[Index(1,0,0)]*(x)*(2-y)*(2-z) +
								 verts[Index(1,0,1)]*(x)*(2-y)*(z) +
								 verts[Index(1,1,0)]*(x)*(y)*(2-z) +
								 verts[Index(1,1,1)]*(x)*(y)*(z) )*.125;

					if (x == 1 || y == 1 || z == 1)
						csg_root->eval(*(vect3d*)&p[x][y][z], p[x][y][z].v[3], g[x][y][z]);
					else
						g[x][y][z] = grad[Index(x>>1,y>>1,z>>1)];
				}
			}
		}

		// create children
		if (guide)
		{
			for (Index i; i < 8; i++)
			{
				for (Index j; j < 8; j++)
				{
					children[i]->verts[j] = p[i.x+j.x][i.y+j.y][i.z+j.z];
					grad[j] = g[i.x+j.x][i.y+j.y][i.z+j.z];
				}

				children[i]->eval(grad, guide->children[i]);
			}
		}
		else
		{
			for (Index i; i < 8; i++)
			{
				children[i] = new TNode();

				for (Index j; j < 8; j++)
				{
					children[i]->verts[j] = p[i.x+j.x][i.y+j.y][i.z+j.z];
					grad[j] = g[i.x+j.x][i.y+j.y][i.z+j.z];
				}

				children[i]->eval(grad, 0);
			}
		}
	}
	else
	{
		// it's a leaf
		//printf("leaf size = %f\n", cellsize);
	}
}


void TNode::defoliate()
{
	for (int i = 0; i < 8; i++)
	{
		delete children[i];
		children[i] = 0;
	}
}

