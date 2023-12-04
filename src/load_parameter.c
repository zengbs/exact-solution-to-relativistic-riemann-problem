#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../include/global.h"

void Load_Parameter() {
#define MAX_STRING_LENGTH 250

  FILE *file = fopen("Input__Parameter", "r");
  if (file == NULL) {
    printf("Error opening file\n");
    exit(1);
  }

  char para_name[MAX_STRING_LENGTH];
  char para_str[MAX_STRING_LENGTH];
  char para_comment[MAX_STRING_LENGTH];
  char line[MAX_STRING_LENGTH];

  while (fgets(line, sizeof(line), file)) {
    //    Empty line or comment
    if (line[0] == '\n' || line[0] == '#') continue;

    //    Parse the line into PARAMETER_NAME, VALUE, and COMMENTS
    int result =
        sscanf(line, "%s %s %[^\n]", para_name, para_str, para_comment);
    if (result < 2) {
      printf("ERROR : Wrong Input__Parameter format!\n");
      continue;
    }

    //    Load the parameter
    if (strcmp(para_name, "L_DENS") == 0) {
      sscanf(para_str, "%lf", &L_DENS);
    } else if (strcmp(para_name, "L_VELX") == 0) {
      sscanf(para_str, "%lf", &L_VELX);
    } else if (strcmp(para_name, "L_PRES") == 0) {
      sscanf(para_str, "%lf", &L_PRES);
    } else if (strcmp(para_name, "R_DENS") == 0) {
      sscanf(para_str, "%lf", &R_DENS);
    } else if (strcmp(para_name, "R_VELX") == 0) {
      sscanf(para_str, "%lf", &R_VELX);
    } else if (strcmp(para_name, "R_PRES") == 0) {
      sscanf(para_str, "%lf", &R_PRES);
    } else if (strcmp(para_name, "DT") == 0) {
      sscanf(para_str, "%lf", &DT);
    } else if (strcmp(para_name, "END_T") == 0) {
      sscanf(para_str, "%lf", &END_T);
    } else if (strcmp(para_name, "L_X") == 0) {
      sscanf(para_str, "%lf", &L_X);
    } else if (strcmp(para_name, "R_X") == 0) {
      sscanf(para_str, "%lf", &R_X);
    } else if (strcmp(para_name, "N_CELL") == 0) {
      sscanf(para_str, "%d", &N_CELL);
    } else {
      printf("Unknown parameter %s !!\n", para_name);
    }  // if ( strcmp(string, "L_Dens") == 0 ) ... else ...
  }

  fclose(file);

#undef MAX_STRING_LENGTH

}  // FUNCTION : Load_Parametrer
