#include <stdio.h>
#include "messages.h"

int f_error(const char *input){
	printf("%s error: %s %s \n",KRED,input,KNRM);
	return 0;
}

int f_message(const char *input){
	printf("%s %s %s \n", KBLU,input,KNRM);
	return 0;
}
int f_usage(){
	printf("%s Using \n  %s",KGRN,KNRM);
	return 0;
}
