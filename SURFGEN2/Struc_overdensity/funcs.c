# include "MarchingCube.h"

double **allocate_double_2d(int N1, int N2)
{
  double **xxa;
  int ii;

  xxa=(double **)malloc(N1 *  sizeof(double*));
  for(ii=0; ii<N1; ii++)
    xxa[ii]=(double *) malloc (N2 * sizeof(double));
  return(xxa);
}

void deallocate_double_2d(double **arr2D, int rows)
{
  int i;

  for(i=0;i<rows;i++)
    {
      free(arr2D[i]);
    }

  free(arr2D);
}

double ABS(double x)
{
  if(x<0)
    return(-1.*x);
  else
    return(x);
}

int size_of_file(char *infile)
{
  int wc;
  FILE *fp;
  char c;

  fp=fopen(infile, "r");

  wc=0;

  while ((c = getc(fp)) != EOF)
    {

      if (c == '\n')
        wc++;
    }

  fclose(fp);
  return(wc);
}

