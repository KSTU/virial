#define MAXATOM 100
#define MAXCHAR 20
#define MAXFILENAME 150
#ifndef SIM_H
#define SIM_H
	typedef struct simulation{
		float final_temp;
		float initial_temp;
		float delta_temp;
		char substance_file_name[MAXFILENAME];
		char substance_top_name[MAXFILENAME];
		char comput[MAXCHAR];
		char type[MAXCHAR];
	}simulation;

	typedef struct molecula{
		char title[MAXCHAR];
		int atom_num;
		int atom_id[MAXATOM];
		float x[MAXATOM];
		float y[MAXATOM];
		float z[MAXATOM];
		float vx[MAXATOM];
		float vy[MAXATOM];
		float vz[MAXATOM];
		char mol_name[MAXCHAR];
		char atom_name[MAXATOM][MAXCHAR];
	}molecula;

	typedef struct topology{
		int number;
		char atom_name[MAXCHAR][MAXATOM];
		float sigma[MAXATOM];
		float epsilon[MAXATOM];
		float q[MAXATOM];
		float mass[MAXATOM];
	}topology;
#endif
