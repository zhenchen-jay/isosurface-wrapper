#pragma once

#include "iso_common.h"
#include "global.h"
#include "array2d.h"
#include "vect.h"
#include "qefnorm.h"
#include "cube_arrays.h"
#include <float.h>
#include "rootfind.h"
#include "csgtree.h"

extern int tree_cells;

struct TNode
{
	TNode()
	{
		children[0] = children[1] = children[2] = children[3] = 
			children[4] = children[5] = children[6] = children[7] = 0;
	}
	~TNode()
	{
		defoliate();
	}

	vect4f verts[8];
	vect4f edges[12];
	vect4f faces[6];
	vect4f node;

	TNode *children[8];

	bool changesSign(vect4f *verts, vect4f *edges, vect4f *faces, vect4f &node);

	void eval(vect3f *grad, TNode *guide);
	bool changesSign(){return changesSign(verts, edges, faces, node);}


	void defoliate();

	bool contains(vect3f &p)
	{
		for (int i = 0; i < 3; i++)
		{
			if (p[i] < verts[0][i] || p[i] > verts[7][i])
				return false;
		}
		return true;
	}

	/*bool is_outside()
	{
		if (verts[7][0] - verts[0][0] < 1.5)
		{
			for (int i = 0; i < 3; i++)
			{
				if (verts[0][i] < 0 || verts[7][i] > 1)
					return true;
			}
		}

		return false;
	}*/
	bool is_outside(){return is_outside((vect3f&)verts[0], (vect3f&)verts[7]);}

	bool is_outside(vect3f &mine, vect3f &maxe)
	{
		/*if (maxe[0] - mine[0] < 1.5)
		{
			for (int i = 0; i < 3; i++)
			{
				if (mine[i] < 0 || maxe[i] > 1)
					return true;
			}
		}*/

		return false;
	}

	bool is_leaf()
	{
		assert(this != 0); 
		return children[0] == 0;
		//if (children[0] == 0)
		//	return true;
		//return is_outside();
	}

	void prune(int depth = 0)
	{
		if (depth >= 2)
		{
			bool all_leaves = true;
			for (int i = 0; i < 8; i++)
			{
				if (children[i] && children[i]->children[0])
				{
					all_leaves = false;
					break;
				}
			}
			if (all_leaves)
				defoliate();
		}

		if (children[0])
			for (int i = 0; i < 8; i++)
				children[i]->prune(depth+1);
	}


	// NODES
	void vertNode(vect4f &p, vect3f *grad, float &qef_error, vect3f *ev = 0, int *ns_size = 0)
	{
		// debug
		vect3f mid_debug = 0;
		for (int i = 0; i < 8; i++)
		{
			vect4f &p4 = verts[i];
			mid_debug += p4;
		}
		mid_debug *= .125;
		bool do_debug = false;

		// build QEF
		//float cellsize = verts[7][0] - verts[0][0];
		QEFNormal<double, 4> q;
		q.zero();

		vect4f mid = 0;

		vector<vect4f> plane_norms, plane_pts;

		for (int x = 0; x <= OVERSAMPLE_QEF; x++)
		{
			for (int y = 0; y <= OVERSAMPLE_QEF; y++)
			{
				for (int z = 0; z <= OVERSAMPLE_QEF; z++)
				{
					if (x == 0 && y == 1 && z == 1)
					{
						int test = 1;
					}

					vect4f p;
					p[0] = (1 - float(x)/OVERSAMPLE_QEF)*verts[0][0] + (float(x)/OVERSAMPLE_QEF)*verts[7][0];
					p[1] = (1 - float(y)/OVERSAMPLE_QEF)*verts[0][1] + (float(y)/OVERSAMPLE_QEF)*verts[7][1];
					p[2] = (1 - float(z)/OVERSAMPLE_QEF)*verts[0][2] + (float(z)/OVERSAMPLE_QEF)*verts[7][2];

					vect5f pl;
					csg_root->eval((vect3f&)p, p[3], (vect3f&)pl);
					pl[3] = -1;
					pl[4] = -(p[0]*pl[0] + p[1]*pl[1] + p[2]*pl[2]) + p[3]; // -p*n

					q.combineSelf(vect5d(pl).v);

					ArrayWrapper<double, 4> tmpA;
					double tmpB[4];

					for (int i = 0; i < 4; i++)
					{
						int index = ((2 * 4 + 3 - i) * i) / 2;
						for (int j = i; j < 4; j++)
						{
							tmpA.data[i][j] = q.data[index + j - i];
							tmpA.data[j][i] = tmpA.data[i][j];
						}

						tmpB[i] = -q.data[index + 4 - i];
					}

					mid += p;

					plane_pts.push_back(p);
					plane_norms.push_back(vect4f(pl[0], pl[1], pl[2], -1));
				}
			}
		}
		mid /= (OVERSAMPLE_QEF+1)*(OVERSAMPLE_QEF+1)*(OVERSAMPLE_QEF+1);

		//** calc minimizer
		if (do_debug)
			int alala=48932;

		// build system to solve
		const int n = 4;
		ArrayWrapper<double, n> A;
		double B [ n ];

		for (int i = 0; i < n; i++ )
		{
			int index = ( ( 2 * n + 3 - i ) * i ) / 2;
			for (int j = i; j < n; j++ )
			{
				A.data [ i ] [ j ] = q.data [ index + j - i ];
				A.data [ j ] [ i ] = A.data [ i ] [ j ];
			}

			B [ i ] = -q.data [ index + n - i ];
		}

		// minimize QEF constrained to cell
		const float border = BORDER * (verts[7][0]-verts[0][0]);
		bool is_out = true;
		double err = 1e30;
		vect3f mine(verts[0][0] + border, verts[0][1] + border, verts[0][2] + border);
		vect3f maxe(verts[7][0] - border, verts[7][1] - border, verts[7][2] - border);

		for (int cell_dim = 3; cell_dim >= 0 && is_out; cell_dim--)
		{
			if (cell_dim == 3)
			{
				// find minimal point
				vect4d rvalue;
				ArrayWrapper<double, n> inv;

				::matInverse<double, n> ( A, inv);

				for (int i = 0; i < n; i++ )
				{
					rvalue [ i ] = 0;
					for (int j = 0; j < n; j++ )
						rvalue [ i ] += inv.data [ j ] [ i ] * B [ j ];
				}

				p(rvalue[0], rvalue[1], rvalue[2], rvalue[3]);

				// check bounds
				if (p[0] >= mine[0] && p[0] <= maxe[0] &&
					p[1] >= mine[1] && p[1] <= maxe[1] &&
					p[2] >= mine[2] && p[2] <= maxe[2])
				{
					is_out = false;
					err = calcError(p, plane_norms, plane_pts);
				}
			}
			else if (cell_dim == 2)
			{
				for (int face = 0; face < 6; face++)
				{
					int dir = face / 2;
					int side = face % 2;
					vect3f corners[2] = {mine, maxe};

					// build constrained system
					ArrayWrapper<double, n+1> AC;
					double BC[n+1];
					for (int i = 0; i < n+1; i++)
					{
						for (int j = 0; j < n+1; j++)
						{
							AC.data[i][j] = (i < n && j < n ? A.data[i][j] : 0);
						}
						BC[i] = (i < n ? B[i] : 0);
					}

					AC.data[n][dir] = AC.data[dir][n] = 1;
					BC[n] = corners[side][dir];

					// find minimal point
					double rvalue[n+1];
					ArrayWrapper<double, n+1> inv;

					::matInverse<double, n+1> ( AC, inv);

					for (int i = 0; i < n+1; i++ )
					{
						rvalue [ i ] = 0;
						for (int j = 0; j < n+1; j++ )
							rvalue [ i ] += inv.data [ j ] [ i ] * BC [ j ];
					}

					vect4f pc(rvalue[0], rvalue[1], rvalue[2], rvalue[3]);
				
					// check bounds
					int dp = (dir+1)%3;
					int dpp = (dir+2)%3;
					if (pc[dp] >= mine[dp] && pc[dp] <= maxe[dp] &&
						pc[dpp] >= mine[dpp] && pc[dpp] <= maxe[dpp])
					{
						is_out = false;
						double e = calcError(pc, plane_norms, plane_pts);
						if (e < err)
						{
							err = e;
							p = pc;
						}
					}
				}
			}
			else if (cell_dim == 1)
			{
				for (int edge = 0; edge < 12; edge++)
				{
					int dir = edge / 4;
					int side = edge % 4;
					vect3f corners[2] = {mine, maxe};

					// build constrained system
					ArrayWrapper<double, n+2> AC;
					double BC[n+2];
					for (int i = 0; i < n+2; i++)
					{
						for (int j = 0; j < n+2; j++)
						{
							AC.data[i][j] = (i < n && j < n ? A.data[i][j] : 0);
						}
						BC[i] = (i < n ? B[i] : 0);
					}

					int dp = (dir+1)%3;
					int dpp = (dir+2)%3;
					AC.data[n][dp] = AC.data[dp][n] = 1;
					AC.data[n+1][dpp] = AC.data[dpp][n+1] = 1;
					BC[n] = corners[side&1][dp];
					BC[n+1] = corners[side>>1][dpp];

					// find minimal point
					double rvalue[n+2];
					ArrayWrapper<double, n+2> inv;

					::matInverse<double, n+2> ( AC, inv);

					for (int i = 0; i < n+2; i++ )
					{
						rvalue [ i ] = 0;
						for (int j = 0; j < n+2; j++ )
							rvalue [ i ] += inv.data [ j ] [ i ] * BC [ j ];
					}

					vect4f pc(rvalue[0], rvalue[1], rvalue[2], rvalue[3]);
				
					// check bounds
					if (pc[dir] >= mine[dir] && pc[dir] <= maxe[dir])
					{
						is_out = false;
						double e = calcError(pc, plane_norms, plane_pts);
						if (e < err)
						{
							err = e;
							p = pc;
						}
					}
				}
			}
			else if (cell_dim == 0)
			{
				for (int vertex = 0; vertex < 8; vertex++)
				{
					vect3f corners[2] = {mine, maxe};

					// build constrained system
					ArrayWrapper<double, n+3> AC;
					double BC[n+3];
					for (int i = 0; i < n+3; i++)
					{
						for (int j = 0; j < n+3; j++)
						{
							AC.data[i][j] = (i < n && j < n ? A.data[i][j] : 0);
						}
						BC[i] = (i < n ? B[i] : 0);
					}

					for (int i = 0; i < 3; i++)
					{
						AC.data[n+i][i] = AC.data[i][n+i] = 1;
						BC[n+i] = corners[(vertex>>i)&1][i];
					}
					// find minimal point
					double rvalue[n+3];
					ArrayWrapper<double, n+3> inv;

					::matInverse<double, n+3> ( AC, inv);

					for (int i = 0; i < n+3; i++ )
					{
						rvalue [ i ] = 0;
						for (int j = 0; j < n+3; j++ )
							rvalue [ i ] += inv.data [ j ] [ i ] * BC [ j ];
					}

					vect4f pc(rvalue[0], rvalue[1], rvalue[2], rvalue[3]);
				
					// check bounds
					double e = calcError(pc, plane_norms, plane_pts);
					if (e < err)
					{
						err = e;
						p = pc;
					}
				}
			}
		}

		// eval
		//float qef_val = p[3];
		fun(p);
		//err = fabs(qef_val - p[3]);

		qef_error += err;
	}

	template <class T, class U>
	static double calcError(T &p, U &plane_norms, U &plane_pts)
	{
		assert(plane_norms.size() > 0);
		double err = 0;
		for (int i = 0; i < plane_norms.size(); i++)
		{
			double c = plane_norms[i]*p - plane_norms[i]*plane_pts[i];
			err += c*c;
		}
		return err;
	}

	// FACES
	void vertFace(vect4f &p, vect3f *grad, float &qef_error, int which_face, vect3f *ev = 0)
	{
		//// debug
		//vect3f mid_debug = 0;
		//for (int i = 0; i < 4; i++)
		//{
		//	vect4f &p4 = verts[cube_face2vert[which_face][i]];
		//	mid_debug += p4;
		//}
		//mid_debug *= .25;
		//bool do_debug = false;
		//if (g.cursor == mid_debug)
		//	do_debug = true;

		// build QEF
		//float cellsize = verts[7][0] - verts[0][0];
		int xi = (cube_face2orient[which_face] + 1) % 3;
		int yi = (cube_face2orient[which_face] + 2) % 3;
		int zi = cube_face2orient[which_face];

		QEFNormal<double, 3> q;
		q.zero();

		vect3f mid = 0;

		vector<vect3f> plane_norms, plane_pts;

		vect4f &c0 = verts[cube_face2vert[which_face][0]], &c1 = verts[cube_face2vert[which_face][2]];
		for (int x = 0; x <= OVERSAMPLE_QEF; x++)
		{
			for (int y = 0; y <= OVERSAMPLE_QEF; y++)
			{
				vect4f p4;
				p4[xi] = (1 - float(x)/OVERSAMPLE_QEF)*c0[xi] + (float(x)/OVERSAMPLE_QEF)*c1[xi];
				p4[yi] = (1 - float(y)/OVERSAMPLE_QEF)*c0[yi] + (float(y)/OVERSAMPLE_QEF)*c1[yi];
				p4[zi] = c0[zi];

				vect3f n4;
				csg_root->eval((vect3f&)p4, p4[3], (vect3f&)n4);
				vect3f n3(n4[xi], n4[yi], -1);
				vect3f p3(p4[xi], p4[yi], p4[3]);
				vect4f pl = n3;
				pl[3] = -(p3*n3);	

			q.combineSelf(vect4d(pl).v);

			plane_pts.push_back(p3);
			plane_norms.push_back(n3);

				mid += p3;
			}
		}
		mid /= (OVERSAMPLE_QEF+1)*(OVERSAMPLE_QEF+1);

		// build system to solve
		const int n = 3;
		ArrayWrapper<double, n> A;
		double B [ n ];

		for (int i = 0; i < n; i++ )
		{
			int index = ( ( 2 * n + 3 - i ) * i ) / 2;
			for (int j = i; j < n; j++ )
			{
				A.data [ i ] [ j ] = q.data [ index + j - i ];
				A.data [ j ] [ i ] = A.data [ i ] [ j ];
			}

			B [ i ] = -q.data [ index + n - i ];
		}

		// minimize QEF constrained to cell
		const float border = BORDER * (verts[7][0]-verts[0][0]);
		bool is_out = true;
		double err = 1e30;
		vect2f mine(verts[cube_face2vert[which_face][0]][xi] + border, verts[cube_face2vert[which_face][0]][yi] + border);
		vect2f maxe(verts[cube_face2vert[which_face][2]][xi] - border, verts[cube_face2vert[which_face][2]][yi] - border);
		vect3f p3;

		for (int cell_dim = 2; cell_dim >= 0 && is_out; cell_dim--)
		{
			if (cell_dim == 2)
			{
				// find minimal point
				vect3d rvalue;
				ArrayWrapper<double, n> inv;

				::matInverse<double, n> ( A, inv);

				for (int i = 0; i < n; i++ )
				{
					rvalue [ i ] = 0;
					for (int j = 0; j < n; j++ )
						rvalue [ i ] += inv.data [ j ] [ i ] * B [ j ];
				}

				p3(rvalue[0], rvalue[1], rvalue[2]);

				// check bounds
				if (p3[0] >= mine[0] && p3[0] <= maxe[0] &&
					p3[1] >= mine[1] && p3[1] <= maxe[1])
				{
					is_out = false;
					err = calcError(p3, plane_norms, plane_pts);
				}
			}
			else if (cell_dim == 1)
			{
				for (int edge = 0; edge < 4; edge++)
				{
					int dir = edge / 2;
					int side = edge % 2;
					vect2f corners[2] = {mine, maxe};

					// build constrained system
					ArrayWrapper<double, n+1> AC;
					double BC[n+1];
					for (int i = 0; i < n+1; i++)
					{
						for (int j = 0; j < n+1; j++)
						{
							AC.data[i][j] = (i < n && j < n ? A.data[i][j] : 0);
						}
						BC[i] = (i < n ? B[i] : 0);
					}

					AC.data[n][dir] = AC.data[dir][n] = 1;
					BC[n] = corners[side][dir];

					// find minimal point
					double rvalue[n+1];
					ArrayWrapper<double, n+1> inv;

					::matInverse<double, n+1> ( AC, inv);

					for (int i = 0; i < n+1; i++ )
					{
						rvalue [ i ] = 0;
						for (int j = 0; j < n+1; j++ )
							rvalue [ i ] += inv.data [ j ] [ i ] * BC [ j ];
					}

					vect4f pc(rvalue[0], rvalue[1], rvalue[2]);
				
					// check bounds
					int dp = (dir+1)%2;
					if (pc[dp] >= mine[dp] && pc[dp] <= maxe[dp])
					{
						is_out = false;
						double e = calcError(pc, plane_norms, plane_pts);
						if (e < err)
						{
							err = e;
							p3 = pc;
						}
					}
				}
			}
			else if (cell_dim == 0)
			{
				for (int vertex = 0; vertex < 4; vertex++)
				{
					vect2f corners[2] = {mine, maxe};

					// build constrained system
					ArrayWrapper<double, n+2> AC;
					double BC[n+2];
					for (int i = 0; i < n+2; i++)
					{
						for (int j = 0; j < n+2; j++)
						{
							AC.data[i][j] = (i < n && j < n ? A.data[i][j] : 0);
						}
						BC[i] = (i < n ? B[i] : 0);
					}

					for (int i = 0; i < 2; i++)
					{
						AC.data[n+i][i] = AC.data[i][n+i] = 1;
						BC[n+i] = corners[(vertex>>i)&1][i];
					}
					// find minimal point
					double rvalue[n+2];
					ArrayWrapper<double, n+2> inv;

					::matInverse<double, n+2> ( AC, inv);

					for (int i = 0; i < n+2; i++ )
					{
						rvalue [ i ] = 0;
						for (int j = 0; j < n+2; j++ )
							rvalue [ i ] += inv.data [ j ] [ i ] * BC [ j ];
					}

					vect3f pc(rvalue[0], rvalue[1], rvalue[2]);
				
					// check bounds
					double e = calcError(pc, plane_norms, plane_pts);
					if (e < err)
					{
						err = e;
						p3 = pc;
					}
				}
			}
		}

		/*if (do_debug)
		{
			printf("ourSoln={%f,%f,%f}\n", p3[0], p3[1], p3[2]);

			printf("NMinimize[{\n");
			for (int i = 0; i < plane_norms.size(); i++)
			{
				if (i!=0)
					printf("+");
				printf("({%f,%f,%f}.{x,y,z}-\n{%f,%f,%f}.{%f,%f,%f})^2\n", 
					plane_norms[i][0], plane_norms[i][1], plane_norms[i][2],
					plane_norms[i][0], plane_norms[i][1], plane_norms[i][2],
					plane_pts[i][0], plane_pts[i][1], plane_pts[i][2]);
			}
			printf(",\n");
			printf("x>=%f&&x<=%f\n", mine[0], maxe[0]);
			printf("&&y>=%f&&y<=%f\n", mine[1], maxe[1]);
			printf("},{x,y,z}]\n");
		}*/

		// unproject back into 4 dimensions
		p = verts[cube_face2vert[which_face][0]];
		p[xi] = p3[0];
		p[yi] = p3[1];
		fun(p);

		qef_error += err;

	}


	// EDGES
	void vertEdge(vect4f &p, vect3f *grad, float &qef_error, int which_edge)
	{
		//// debug
		//vect3f mid_debug = 0;
		//for (int i = 0; i < 2; i++)
		//{
		//	vect4f &p4 = verts[cube_edge2vert[which_edge][i]];
		//	mid_debug += p4;
		//}
		//mid_debug *= .5;
		//bool do_debug = false;
		//if (g.cursor == mid_debug)
		//	do_debug = true;

		// calc QEF
		//float cellsize = verts[7][0] - verts[0][0];
		int xi = cube_edge2orient[which_edge];
		int yi = (cube_edge2orient[which_edge] + 1) % 3;
		int zi = (cube_edge2orient[which_edge] + 2) % 3;

		QEFNormal<double, 2> q;
		q.zero();

		vect2f mid = 0;

		vector<vect2f> plane_norms, plane_pts;

		vect4f &c0 = verts[cube_edge2vert[which_edge][0]], &c1 = verts[cube_edge2vert[which_edge][1]];
		for (int i = 0; i <= OVERSAMPLE_QEF; i++)
		{
			vect4f p4;
			p4[xi] = (1 - float(i)/OVERSAMPLE_QEF)*c0[xi] + (float(i)/OVERSAMPLE_QEF)*c1[xi];
			p4[yi] = c0[yi];
			p4[zi] = c0[zi];

			vect3f g3;
			csg_root->eval((vect3f&)p4, p4[3], (vect3f&)g3);
			vect2f n2(g3[xi], -1);
			vect2f p2(p4[xi], p4[3]);
			vect3f pl = n2;
			pl[2] = -(p2*n2);

			q.combineSelf(vect3d(pl).v);

			plane_norms.push_back(n2);
			plane_pts.push_back(p2);

			mid += p2;
		}
		mid /= OVERSAMPLE_QEF+1;

		// build system to solve
		const int n = 2;
		ArrayWrapper<double, n> A;
		double B [ n ];

		for (int i = 0; i < n; i++ )
		{
			int index = ( ( 2 * n + 3 - i ) * i ) / 2;
			for (int j = i; j < n; j++ )
			{
				A.data [ i ] [ j ] = q.data [ index + j - i ];
				A.data [ j ] [ i ] = A.data [ i ] [ j ];
			}

			B [ i ] = -q.data [ index + n - i ];
		}

		// minimize QEF constrained to cell
		const float border = BORDER * (verts[7][0]-verts[0][0]);
		bool is_out = true;
		double err = 1e30;
		const float vmin = verts[cube_edge2vert[which_edge][0]][xi] + border;
		const float vmax = verts[cube_edge2vert[which_edge][1]][xi] - border;
		vect2f p2;

		for (int cell_dim = 1; cell_dim >= 0 && is_out; cell_dim--)
		{
			if (cell_dim == 1)
			{
				// find minimal point
				vect3d rvalue;
				ArrayWrapper<double, n> inv;

				::matInverse<double, n> ( A, inv);

				for (int i = 0; i < n; i++ )
				{
					rvalue [ i ] = 0;
					for (int j = 0; j < n; j++ )
						rvalue [ i ] += inv.data [ j ] [ i ] * B [ j ];
				}

				p2(rvalue[0], rvalue[1]);

				// check bounds
				if (p2[0] >= vmin && p2[0] <= vmax)
				{
					is_out = false;
					err = calcError(p2, plane_norms, plane_pts);
				}
			}
			else if (cell_dim == 0)
			{
				for (int vertex = 0; vertex < 2; vertex++)
				{
					float corners[2] = {vmin, vmax};

					// build constrained system
					ArrayWrapper<double, n+1> AC;
					double BC[n+1];
					for (int i = 0; i < n+1; i++)
					{
						for (int j = 0; j < n+1; j++)
						{
							AC.data[i][j] = (i < n && j < n ? A.data[i][j] : 0);
						}
						BC[i] = (i < n ? B[i] : 0);
					}

					for (int i = 0; i < 1; i++)
					{
						AC.data[n+i][i] = AC.data[i][n+i] = 1;
						BC[n+i] = corners[vertex>>i];
					}
					// find minimal point
					double rvalue[n+1];
					ArrayWrapper<double, n+1> inv;

					::matInverse<double, n+1> ( AC, inv);

					for (int i = 0; i < n+1; i++ )
					{
						rvalue [ i ] = 0;
						for (int j = 0; j < n+1; j++ )
							rvalue [ i ] += inv.data [ j ] [ i ] * BC [ j ];
					}

					vect2f pc(rvalue[0], rvalue[1]);
				
					// check bounds
					double e = calcError(pc, plane_norms, plane_pts);
					if (e < err)
					{
						err = e;
						p2 = pc;
					}
				}
			}
		}

		
		/*if (do_debug)
		{
			printf("ourSoln={%f,%f}\n", p2[0], p2[1]);

			printf("NMinimize[{\n");
			for (int i = 0; i < plane_norms.size(); i++)
			{
				if (i!=0)
					printf("+");
				printf("({%f,%f}.{x,y}-\n{%f,%f}.{%f,%f})^2\n", 
					plane_norms[i][0], plane_norms[i][1],
					plane_norms[i][0], plane_norms[i][1],
					plane_pts[i][0], plane_pts[i][1]);
			}
			printf(",\n");
			printf("x>=%f&&x<=%f\n", vmin, vmax);
			printf("},{x,y}]\n");
		}*/

		// unproject back into 4 dimensions
		p = verts[cube_edge2vert[which_edge][0]];
		p[xi] = p2[0];
		fun(p);

		qef_error += err;
		
	}
};

void gen_iso_ours();
