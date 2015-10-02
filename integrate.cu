#include <stdio.h>
#include "input.h"
#include "messages.h"
#include "simulation.h"
#include "integrate.h"

int prop_boundary(molecula *sub,topology *top){
	FILE *temp_file;
	float r_cur;
	float r_min;
	float r_delta;
	int dirrection; //1 -- up, -1 -- down
	float mix_sigma;
	float mix_epsilon;
	float mix_q;
	int iter;
	float f1,f2;

	r_cur=abs(sub->x[1]-sub->x[0]);
	dirrection=1;
	r_delta=r_cur/1000;
	
	//get potential minimum
	mix_epsilon=mix_geom(top->sigma[sub->atom_id[0]],top->sigma[sub->atom_id[1]]);
	mix_sigma=mix_ariph(top->epsilon[sub->atom_id[0]],top->epsilon[sub->atom_id[1]]);
	mix_q=mix_charge(top->q[sub->atom_id[0]],top->q[sub->atom_id[1]]);
	iter=0;
	while(iter<100){
		f1=potential(mix_sigma,mix_epsilon,mix_q,r_cur);
		f2=potential(mix_sigma,mix_epsilon,mix_q,r_cur+dirrection*r_delta);
		if(f1>f2){
			dirrection=-dirrection;
			r_delta=r_delta/2.0;
		}
		else{
			r_cur=r_cur+dirrection*r_delta;
		}
	}
	//get probability distribution
	
	
	//write to file
	temp_file=fopen("prop.out","w");
	
	fclose(temp_file);
	return 0;
}

float potential(float sig,float eps, float q, float r){
	float p;
	float sr;
	const float k=1.08;	//
	
	sr=sig/r;
	sr=sr*sr;	//2
	sr=sr*sr*sr;	//6
	p=4.0*eps*(sr*sr-sr)+k*q/r;
	
	return 0.2;
}

float mix_ariph(float s1,float s2){
	return (s1+s2)/2.0;
}

float mix_geom(float s1, float s2){
	return sqrt(s1*s2);
}
float mix_charge(float q1,float q2){
	return q1*q2;
}
