#include <stdio.h>
#include "messages.h"
#include "input.h"
#include "simulation.h"


int main(int argc, char *argv[]){
	int DeviceCount;
	cudaDeviceProp dp;
	simulation sp;
	molecula substance;
	topology top;
	int i;
	//
	
	//
	cudaGetDeviceCount(&DeviceCount);
	printf("Found %d device \n",DeviceCount);
	for (int device =0; device<DeviceCount;device++){
		cudaGetDeviceProperties(&dp,device);
		printf("Clock rate : %d \n", dp.clockRate);
		printf("Max thread dimention %d %d %d \n", dp.maxThreadsDim[0],dp.maxThreadsDim[1],dp.maxThreadsDim[2]);
	}
	if(argc<3){
		f_usage();
		return 0;
	}
	if(check_flag(argc,argv,&sp)!=0){
		f_error("checking program parameters");
		return 1;
	}
	//read initial gro file
	if(read_gro(sp.substance_file_name,&substance)!=0){
		f_error("reading gro file");
		return 1;
	}
	//read topology
	if(read_top(sp.substance_top_name,&top)!=0){
		f_error("reading topology file");
		return 1;
	}
	//get atoms ID
	if(get_atom_id(&substance,&top)!=0){
		f_error("getting atom ID");
		return 1;
	}
	if(strcmp(sp.type,"IL")==0){
		f_message("compute for ionic liquid type");
		if (prop_boundary(&sub,&top)!=0){
			f_error("boundary fail");
			return 1;
		}
		
	}
}
