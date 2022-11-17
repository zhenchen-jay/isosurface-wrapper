#pragma once

#include "vect.h"
#include "global.h"

#include <map>
#include <stdio.h>
#include <Eigen/Dense>
#include <igl/AABB.h>

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

struct GeneralShape : public CSGNode
{
	GeneralShape(Eigen::MatrixXf& V, Eigen::MatrixXi& F, float eps)
	{
		_V = V;
		_F = F;

		// shrink and translate mesh to the box [0, 0, 0], and [1, 1, 1]
		Eigen::Vector3f minCorner = V.colwise().minCoeff();
		Eigen::Vector3f maxCorner = V.colwise().maxCoeff();
		_center = (maxCorner + minCorner) / 2;
		Eigen::Vector3f tarCenter;
		tarCenter << 0.5, 0.5, 0.5;

		float diag = (maxCorner - minCorner).norm();
		float max = (maxCorner - minCorner).maxCoeff();
		_ratio = float(1.) / (max + 4 * eps * diag);

		for (int i = 0; i < _V.rows(); i++)
			_V.row(i) = _V.row(i) - _center.transpose();		// move center to [0, 0, 0]
		_V = _V * _ratio;
		_eps = eps * diag * _ratio;		// within [-1/2, -1/2, -1/2] and [1/2, 1/2, 1/2]
		for (int i = 0; i < _V.rows(); i++)
			_V.row(i) = _V.row(i) + tarCenter.transpose();	// move center to [1/2, 1/2, 1/2]

		minCorner = _V.colwise().minCoeff();
		maxCorner = _V.colwise().maxCoeff();

		std::cout << minCorner.transpose() << " " << maxCorner.transpose() << ", eps: " << _eps << ", input: " << eps * diag << std::endl;

		_tree.init(_V, _F);

		float lower_bound = std::numeric_limits<float>::min();
		float upper_bound = std::numeric_limits<float>::max();
		const float max_abs = std::max(std::abs(lower_bound), std::abs(upper_bound));
		up_sqr_d = std::pow(max_abs, 2.0);
		low_sqr_d = std::pow(std::max(max_abs - (upper_bound - lower_bound), (float)0.0), 2.0);
	}

	Eigen::MatrixXf _V;
	Eigen::MatrixXi _F;
	Eigen::Vector3f _center;
	float _ratio;
	igl::AABB<Eigen::MatrixXf, 3> _tree;
	float _eps;
	float up_sqr_d;
	float low_sqr_d;

	virtual void eval(vect3f& pt, float& val, vect3f& grad)
	{
		Eigen::Matrix<float, 1, 3> p, c;
		p << pt[0], pt[1], pt[2];
		int i;
		float sqrd = _tree.squared_distance(_V, _F, p, low_sqr_d, up_sqr_d, i, c);

		if (sqrd >= up_sqr_d || sqrd < low_sqr_d)
			val = std::numeric_limits<float>::quiet_NaN();
		else
			val = std::sqrt(sqrd) - _eps;

		grad[0] = p[0] - c[0];
		grad[1] = p[1] - c[1];
		grad[2] = p[2] - c[2];

		float l = grad.length();
		if (l > 1e-9)
			grad = grad / l;
		else
			grad(0, 0, 1);
	}
};
