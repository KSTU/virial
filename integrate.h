#ifndef INTE_H
#define INTE_H
	#include "simulation.h"
	int prop_boundary(molecula *sub,topology *top);
	float potential(float sig,float eps, float q, float r);
	float mix_ariph(float s1,float s2);
	float mix_geom(float s1, float s2);
	float mix_charge(float q1,float q2);
#endif
