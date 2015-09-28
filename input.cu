#include <stdio.h>
#include "input.h"
#include "messages.h"
#include "simulation.h"
//checking flags
int check_flag(int flag_num, char **param,simulation *sp){
	int i;
	FILE *input_file;
	char temp_string[MAXCHAR];
	float temp_float;
	//
	for(i=0;i<flag_num;i++){
		if(param[i][0]=='-'){
			//simulation parameters
			if(param[i][1]=='p'){
				if(param[i+1][0]!='-'){
					input_file=fopen(param[i+1],"r");
					if(input_file!=NULL){
						f_message("reading parameters");

						fscanf(input_file,"%f",&temp_float);
						sp->initial_temp=temp_float;
						printf(" Initial temperature: %f \n",sp->initial_temp);

						fscanf(input_file,"%f",&temp_float);
						sp->final_temp=temp_float;
						printf(" Final temperature: %f \n",sp->final_temp);

						fscanf(input_file,"%f",&temp_float);
						sp->delta_temp=temp_float;
						printf(" Temperature step: %f \n",sp->delta_temp);
						
						fscanf(input_file,"%s",sp->substance_file_name);
						printf(" File name: %s \n", sp->substance_file_name);

						fscanf(input_file,"%s",sp->substance_top_name);
						printf(" Topology file name: %s \n", sp->substance_top_name);

						fscanf(input_file,"%s",sp->comput);
						printf(" Compute on: %s \n", sp->comput);
						
						fscanf(input_file,"%s",sp->type);
						printf(" Type: %s \n", sp->type);
					}
					else{
						f_error("open parameters file");
					}
					fclose(input_file);
				}
				else{
					f_error("missing simulation parameters file");
				}
			}
			else{
				temp_string[0]='u';	//crutch
				strcat(temp_string, "nknown parameter ");
				f_error(strcat(temp_string,param[i]));
			}
		}
	}
	return 0;
}
int read_gro(const char *input,molecula *sub){
	FILE *temp_file;
	char temp_string[2*MAXFILENAME];
	int temp_int;
	int temp_int2;
	int i;
	
	temp_string[0]='\0';
	strcat(temp_string,"reading gro file: ");
	strcat(temp_string,input);
	f_message(temp_string);

	temp_file=fopen(input,"r");
	if(temp_file==NULL){
		f_error("error in open file");
		return 1;
	}
	fscanf(temp_file,"%20s",sub->title);
	fscanf(temp_file,"%d",&sub->atom_num);
	printf(" atom num %d \n", sub->atom_num);
	for(i=0;i<sub->atom_num;i++){
		fscanf(temp_file,"%5d%5s%5s%5d%f%f%f",&temp_int,sub->mol_name,sub->atom_name[i],&temp_int2,&sub->x[i],&sub->y[i],&sub->z[i]); 
	}
	printf("molecula: %s \n", sub->mol_name);
	for(i=0;i<sub->atom_num;i++){
		printf(" atom type: %s x %f y %f z %f \n",sub->atom_name[i],sub->x[i],sub->y[i],sub->z[i]);
	}
	fclose(temp_file);
	return 0;
}

int line_count(char *input_file_name){
	char c;
	int lines = 0;
	FILE *checking_file;
	//
	checking_file=fopen(input_file_name,"r");
	while ((c=fgetc(checking_file))!=EOF){
		if(c=='\n') ++lines;
	}
	fclose(checking_file);
	return lines;
}

int read_top(char *file_name,topology *in){
	FILE *temp_file;
	int i;
	char temp_string[MAXFILENAME];

	in->number=line_count(file_name)-1;
	temp_file=fopen(file_name,"r");
	fscanf(temp_file,"%s",temp_string);
	for (i=0;i<in->number;i++){
		fscanf(temp_file,"%s%f%f%f%f",in->atom_name[i],&in->sigma[i],&in->epsilon[i],&in->q[i],&in->mass[i]);
	}
	fclose(temp_file);
	return 0;
}

int get_atom_id(molecula *sub,topology *top){
	int i,j;

	for(i=0;i<sub->atom_num;i++){
		for(j=0;j<top->number;j++){
			if(strcmp(sub->atom_name[i],top->atom_name[j])==0){
				sub->atom_id[i]=j;
			}
		}
	}
	return 0;
}
