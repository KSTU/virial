// вывод информации на экран
#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KBLU  "\x1B[34m"

#ifndef ERR_H
#define ERR_H
	int f_error(const char *input);
	int f_message(const char *input);
	int f_usage();
#endif 

