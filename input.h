#ifndef INP_H
#define INP_H
	#include "simulation.h"
	int read_gro(const char *input,molecula *sub);
	int check_flag(int flag_num, char **param,simulation *sp);
	int line_count(char *input_file_name);
	int read_top(char *file_name,topology *in);
	int get_atom_id(molecula *sub,topology *top);
#endif
