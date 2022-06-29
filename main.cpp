#include <iostream>
#include <Eigen/Dense>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>

#include "global.h"
#include "iso_common.h"
#include "iso_method_ours.h"
#include "mctable.h"
#include "csgtree.h"

// for visualization
#include "polyscope/polyscope.h"
#include "polyscope/messages.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"


#include <set>

int OVERSAMPLE_QEF = 4;
double BORDER = (1.0 / 16.0);
int DEPTH_MAX = 2; // 7
int DEPTH_MIN = 1; // 4
int FIND_ROOT_DEPTH = 0;


CSGNode *csg_root = 0;
using namespace std;

void writeFile(Mesh &m, string fn)
{
	FILE *f = fopen(fn.c_str(), "w");

#ifdef JOIN_VERTS
	map<TopoEdge, int> vt;

	int vert_num = 0;
	for (int i = 0; i < m.tris.size(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			map<TopoEdge, int>::iterator it = vt.find(m.topoTris[i][j]);
			if (it == vt.end())
			{
				vert_num++;

				vt[m.topoTris[i][j]] = vert_num;
				
				vect3d pos = m.tris[i][j];
				fprintf(f, "v %f %f %f\n", pos[0], pos[1], pos[2]);
			}
		}
	}

	set<vect<2, TopoEdge> > edges; // only for genus
	for (int i = 0; i < m.tris.size(); i++)
	{
		fprintf(f, "f %d %d %d\n", vt[m.topoTris[i][0]], vt[m.topoTris[i][1]], vt[m.topoTris[i][2]]);
		
		for (int j = 0; j < 3; j++)
		{
			vect<2, TopoEdge> e;
			e[0] = m.topoTris[i][j];
			e[1] = m.topoTris[i][(j+1) % 3];
			if (e[1] < e[0])
				swap(e[0], e[1]);
			edges.insert(e);
		}
	}

	// check genus
	int edge_num = edges.size();
	int face_num = m.tris.size();
	printf("verts = %d, edges = %d, faces = %d, genus = %d\n", vert_num, edge_num, face_num, vert_num - edge_num + face_num);
#else
	for (int j = 0; j < m.tris.size(); j++)
	{
		for (int i = 0; i < 3; i++)
			fprintf(f, "v %f %f %f\n", m.tris[j][i][0], m.tris[j][i][1], m.tris[j][i][2]);

		fprintf(f, "f %d %d %d\n", j*3+1, j*3+2, j*3+3);
	}
#endif
	fclose(f);
}

static void convert2EigenMesh(Mesh &m, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
	map<TopoEdge, int> vt;
	std::vector<Eigen::Vector3d> posList = {};
	std::vector<Eigen::Vector3i> faceList = {};

	int vert_num = 0;
	for (int i = 0; i < m.tris.size(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			map<TopoEdge, int>::iterator it = vt.find(m.topoTris[i][j]);
			if (it == vt.end())
			{
				vert_num++;

				vt[m.topoTris[i][j]] = vert_num;

				vect3d pos = m.tris[i][j];
				posList.push_back(Eigen::Vector3d(pos[0], pos[1], pos[2]));
			}
		}
	}

	set<vect<2, TopoEdge> > edges; // only for genus
	for (int i = 0; i < m.tris.size(); i++)
	{
		faceList.push_back(Eigen::Vector3i(vt[m.topoTris[i][0]] - 1, vt[m.topoTris[i][1]] - 1, vt[m.topoTris[i][2]] - 1));

		for (int j = 0; j < 3; j++)
		{
			vect<2, TopoEdge> e;
			e[0] = m.topoTris[i][j];
			e[1] = m.topoTris[i][(j+1) % 3];
			if (e[1] < e[0])
				swap(e[0], e[1]);
			edges.insert(e);
		}
	}

	// check genus
	int edge_num = edges.size();
	int face_num = m.tris.size();
	printf("verts = %d, edges = %d, faces = %d, genus = %d\n", vert_num, edge_num, face_num, vert_num - edge_num + face_num);

	V.resize(posList.size(), 3);
	F.resize(faceList.size(), 3);

	for(int i = 0; i < V.rows(); i++)
		V.row(i) = posList[i];
	for(int i = 0; i < F.rows(); i++)
		F.row(i) = faceList[i];
}

void init()
{
	// define/load function
	if (!csg_root)
	{
		CSGNode *n;

//		n = new CSGMax(new CSGPlane(vect3d(.7,.7,.7 - .05), -~vect3d(1,1,1)),new CSGPlane(vect3d(.7,.7,.7), ~vect3d(1,1,1))); // used for dc topology figure depth 5 inverted
		//n = new CSGMin(new CSGPlane(vect3d(.7,.7,.7), -~vect3d(1,1,1)),new CSGPlane(vect3d(.7,.7,.7 - .05), ~vect3d(1,1,1))); // used for dc topology figure depth 5
//		n = new CSGSphere(vect3d(.5, .5, .5), .4);
//		n = new CSGTorus(vect3d(.5, .5, .5), .3, .1);
//		n = new CSGMin(new CSGSphere(vect3d(.5, .6, .35 - .01), .25),new CSGSphere(vect3d(.5, .5, .65 - .01), .25));
		//n = new CSGMax(new CSGSphere(vect3d(.5, .5, .3 - .05), .2),new CSGSphere(vect3d(.5, .5, .7 - .05), .2));
		//n = new CSGMin(new CSGSphere(vect3d(.5, .5, .65), .25),new CSGNeg(new CSGCylinder(vect3d(.5, .5, .35), vect3d(0,0,1), .15)));
		//n = new CSGMin(	new CSGSphere(vect3d(.5, .5, .5), .4),	new CSGPlane(vect3d(.65, .65, .55), vect3d(0,0,1))	);
		//n = new CSGPlane(vect3d(.65, .65, .65), vect3d(.2,.4,-1));

		// figure of interesting shapes in one cell
//		n = new CSGSphere(vect3d(.85, .25, .25), .1);
//		n = new CSGCylinder(vect3d(.85, .25, .25), vect3d(0,0,1), .1);
//		n = new CSGCylinder(vect3d(.85, .25, .25), ~vect3d(0.1,.7,.7), .1);
//		n = new CSGMax(new CSGPlane(vect3d(.75, .65, .15), -~vect3d(0,.5,1)),new CSGPlane(vect3d(.75, .65, .15), ~vect3d(0,0,1)));
//		n = new CSGMax(
//			new CSGPlane(vect3d(.75, .65, .15), -~vect3d(0,0,1)),
//			new CSGPlane(vect3d(.75, .65, .15), ~vect3d(0.5,0,1))
//			);
//		n = new CSGMax(new CSGMax(
//			new CSGPlane(vect3d(.35, .25, .35), -~vect3d(0,.5,1)),
//			new CSGPlane(vect3d(.35, .25, .35), ~vect3d(0,0,1))),
//			new CSGPlane(vect3d(.35, .25, .35), -~vect3d(.5,0,1)) );
		
		// draw box with things cut out (shifted) (mechanical part)

//		vect3d shift(.062151346, .0725234, .0412);
//		CSGNode *box =
//					new CSGMin(
//						new CSGMin(
//							new CSGMin(
//								new CSGPlane(vect3d(.3, .3, .3)+shift, vect3d(1,0,0)),
//								new CSGPlane(vect3d(.3, .3, .3)+shift, vect3d(0,1,0)) ),
//							new CSGMin(
//								new CSGPlane(vect3d(.3, .3, .3)+shift, vect3d(0,0,1)),
//								new CSGPlane(vect3d(.7, .7, .7)+shift, vect3d(-1,0,0)) )
//								),
//						new CSGMin(
//							new CSGPlane(vect3d(.7, .7, .7)+shift, vect3d(0,-1,0)),
//							new CSGPlane(vect3d(.7, .7, .7)+shift, vect3d(0,0,-1)) )
//					);
//		n = new CSGMin(new CSGNeg(new CSGCylinder(vect3d(.5, .5, .5)+shift, vect3d(1,0,0), .15)),
//					new CSGMin(new CSGNeg(new CSGCylinder(vect3d(.5, .5, .5)+shift, vect3d(0,1,0), .15)),
//						new CSGMax(box, new CSGCylinder(vect3d(.5, .5, .5)+shift, vect3d(0,0,1), .15) ) )) ;
//		n = new CSGMin(n, new CSGMin(new CSGPlane(vect3d(0, 0, .9)+shift, vect3d(0,0,-1)), new CSGPlane(vect3d(0, 0, .1)+shift, vect3d(0,0,1))));
		n = new CSGMin(
			new CSGMin(
				new CSGMin(
					new CSGPlane(vect3d(.3, .3, .3), vect3d(1, 0, 0)),
					new CSGPlane(vect3d(.3, .3, .3), vect3d(0, 1, 0))),
				new CSGMin(
					new CSGPlane(vect3d(.3, .3, .3), vect3d(0, 0, 1)),
					new CSGPlane(vect3d(.7, .7, .7) , vect3d(-1, 0, 0)))
			),
			new CSGMin(
				new CSGPlane(vect3d(.7, .7, .7), vect3d(0, -1, 0)),
				new CSGPlane(vect3d(.7, .7, .7), vect3d(0, 0, -1)))
		);

		// $$$$$$ add it $$$$$$
		csg_root = n;
	}

	// generate surfaces
//	initMCTable();

	gen_iso_ours(); // pregenerate tree
	//gen_iso_ours();

//	writeFile(g.ourMesh, "output.obj");
}

int main(int argc, char **argv)
{
	::init();
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	convert2EigenMesh(g.ourMesh, V, F);
	igl::writeOBJ("output.obj", V, F);

	// Options
	//polyscope::options::autocenterStructures = false;
	polyscope::view::windowWidth = 1024;
	polyscope::view::windowHeight = 1024;

	// Initialize polyscope
	polyscope::init();
	// Register the mesh with Polyscope
	polyscope::registerSurfaceMesh("iso-surface", V, F);

	// Show the gui
	polyscope::show();

	return 0;
}
