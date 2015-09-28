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
	mix_epsilon=0.1;
	mix_sigma=0.1;
	mix_q=0.1;
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
	
	fopen(temp_file);
	return 0;
}
