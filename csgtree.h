#pragma once

#include "vect.h"
#include "global.h"

#include <map>
#include <stdio.h>

using namespace std;

struct CSGNode
{
	virtual void eval(vect3f &pt, float &val, vect3f &grad) = 0;
};

struct CSGMax : public CSGNode
{
	CSGMax(CSGNode *a, CSGNode *b)
	{
		n1 = a;
		n2 = b;
	}

	virtual void eval(vect3f &pt, float &val, vect3f &grad)
	{
		vect3f g1, g2;
		float v1, v2;

		n1->eval(pt, v1, g1);
		n2->eval(pt, v2, g2);

		if (v1 > v2)
		{
			val = v1;
			grad = g1;
		}
		else
		{
			val = v2;
			grad = g2;
		}
	}

	CSGNode *n1, *n2;
};

struct CSGMin : public CSGNode
{
	CSGMin(CSGNode *a, CSGNode *b)
	{
		n1 = a;
		n2 = b;
	}

	virtual void eval(vect3f &pt, float &val, vect3f &grad)
	{
		vect3f g1, g2;
		float v1, v2;

		n1->eval(pt, v1, g1);
		n2->eval(pt, v2, g2);

		if (v1 < v2)
		{
			val = v1;
			grad = g1;
		}
		else
		{
			val = v2;
			grad = g2;
		}
	}

	CSGNode *n1, *n2;
};

struct CSGDiff : public CSGNode
{
	CSGDiff(CSGNode *a, CSGNode *b)
	{
		n1 = a;
		n2 = b;
	}

	virtual void eval(vect3f &pt, float &val, vect3f &grad)
	{
		vect3f g1, g2;
		float v1, v2;

		n1->eval(pt, v1, g1);
		n2->eval(pt, v2, g2);
		v2 *= -1;
		g2 *= -1;

		if (v1 < v2)
		{
			val = v1;
			grad = g1;
		}
		else
		{
			val = v2;
			grad = g2;
		}
	}

	CSGNode *n1, *n2;
};

struct CSGNeg : public CSGNode
{
	CSGNeg(CSGNode *a)
	{
		n = a;
	}

	virtual void eval(vect3f &pt, float &val, vect3f &grad)
	{
		n->eval(pt, val, grad);
		val = -val;
		grad = -grad;
	}

	CSGNode *n;
};

struct CSGPlane : public CSGNode
{
	CSGPlane(vect3f p, vect3f n)
	{
		pos = p;
		norm = ~n;
	}

	vect3f pos, norm;

	virtual void eval(vect3f &pt, float &val, vect3f &grad)
	{
		grad = norm;
		val = (pt - pos) * norm;
	}
};

struct CSGSphere : public CSGNode
{
	CSGSphere(vect3f p, float r)
	{
		pos = p;
		radius = r;
	}

	vect3f pos;
	float radius;

	virtual void eval(vect3f &pt, float &val, vect3f &grad)
	{
		vect3f d = pt - pos;
		float len = d.length();
		if (len > 1e-9)
			grad = -d / len;
		else
			grad(0,0,1);
		val = radius - len;
	}
};

struct CSGTorus : public CSGNode
{
	CSGTorus(vect3f p, float r1, float r2)
	{
		pos = p;
		R1 = r1*r1;
		R2 = r2*r2;
	}

	vect3f pos;
	float R1, R2;

	virtual void eval(vect3f &pt, float &val, vect3f &grad)
	{
		float x = pt[0] - pos[0];
		float y = pt[1] - pos[1];
		float z = pt[2] - pos[2];

		float a = x*x + y*y + z*z + R1 - R2;
		val = -a*a + 4*R1*(x*x + y*y);

		grad[0] = 8*R1*x - 4*x*a;
		grad[1] = 8*R1*y - 4*y*a;
		grad[2] = -4*z*a;
	}
};

struct CSGCylinder : public CSGNode
{
	CSGCylinder(vect3f p, vect3f d, float r)
	{
		dir = ~d;
		pos = p;
		radius = r;
	}

	vect3f pos, dir;
	float radius;

	virtual void eval(vect3f &pt, float &val, vect3f &grad)
	{
		vect3f v = pt - pos;
		vect3f d = v - dir*(v*dir);
		float l = d.length();
		if (l > 1e-9)
			grad = -d / l;
		else
			grad(0,0,1);
		val = radius - l;
	}
};
