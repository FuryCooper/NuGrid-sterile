#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "vars.h"

void read_parameterfile(char *fname)
{
#define FLOAT 1
#define STRING 2
#define INT 3
#define MAXTAGS 300

  FILE *fd;
  char buf[200], buf1[200], buf2[200], buf3[200];
  int i, j, nt;
  int id[MAXTAGS];
  void *addr[MAXTAGS];
  char *ret, tag[MAXTAGS][50];
  int errorFlag = 0;

  ThisTask = 0;

  /* read parameter file on all processes for simplicty */

  nt = 0;

  strcpy(tag[nt], "output");
  addr[nt] = output;
  id[nt++] = STRING;
    
    strcpy(tag[nt], "input_path");
    addr[nt] = input_path;
    id[nt++] = STRING;

  strcpy(tag[nt], "snapno");
  addr[nt] = &snapno;
  id[nt++] = INT;
    
    strcpy(tag[nt], "nu_txt");
    addr[nt] = nu_txt;
    id[nt++] = STRING;

    
    strcpy(tag[nt], "fnu");
    addr[nt] = &fnu;
    id[nt++] = FLOAT;
    

  if((fd = fopen(fname, "r")))
    {
      while(!feof(fd))
  {
    buf[0] = 0;
    ret = fgets(buf, 200, fd);

    if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
      continue;

    if(buf1[0] == '%')
      continue;

    for(i = 0, j = -1; i < nt; i++)
      if(strcmp(buf1, tag[i]) == 0)
        {
    j = i;
    tag[i][0] = 0;
    break;
        }

    if(j >= 0)
      {
        switch (id[j])
    {
    case FLOAT:
      *((double *) addr[j]) = atof(buf2);
      break;
    case STRING:
      strcpy(addr[j], buf2);
      break;
    case INT:
      *((int *) addr[j]) = atoi(buf2);
      break;
    }
      }
    else
      {
        if(ThisTask == 0)
    fprintf(stdout, "Error in file %s:   Tag '%s' not allowed or multiple defined.\n", fname,
      buf1);
        errorFlag = 1;
      }
  }
      fclose(fd);

    }
  else
    {
      if(ThisTask == 0)
  fprintf(stdout, "Parameter file %s not found.\n", fname);
      errorFlag = 1;
    }


  for(i = 0; i < nt; i++)
    {
      if(*tag[i])
  {
    if(ThisTask == 0)
      fprintf(stdout, "Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
    errorFlag = 1;
  }
    }

  if(errorFlag)
    {
      MPI_Finalize();
      exit(0);
    }


#undef FLOAT
#undef STRING
#undef INT
#undef MAXTAGS
}
