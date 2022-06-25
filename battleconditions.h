#pragma once

#include <vector>
#include "vect.h"

using namespace std;
struct GeoEntity;
struct GeoTileRef;

struct BattleConditions
{
	int aggressor;
	GeoTileRef *tile;
	vector<GeoEntity*> units[2];

	void fight();
};