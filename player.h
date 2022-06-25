#pragma once

#include "vect.h"

extern float pi;

struct Player
{
	Player()
	{
		pos(0,0,.5);
		pitch = yaw = 0;
	}

	vect3f pos;
	float pitch, yaw;
	vect3f n1, n2, n3, n4;

	vect3f getForward()
	{
		return vect3f(cos(yaw)*cos(pitch), sin(yaw)*cos(pitch), sin(pitch));
	}
	vect3f getRight()
	{
		return vect3f(sin(yaw), -cos(yaw), 0);
	}
	vect3f getUp()
	{
		return getRight()%getForward();
	}
};
