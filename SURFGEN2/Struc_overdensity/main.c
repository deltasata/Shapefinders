/*
 * main.c
 *
 *  Created on: 13-Feb-2014
 *      Author: prakash
 */

# include "MarchingCube.h"

int main(int argc, char **argv)  {
  //double area, vol, imc, genus;
  int NC;// ii, index;
  //double Shape1, Shape2, Shape3, Planarity, Filamentarity;
 // double s1, s2, s3, temp;
  //float INDEX;
  int unused __attribute__((unused));
  //FILE *outfile;
  //int NT=0;
  
  if(argc!=7) {
    printf("usage:\n\t %s <inputfilename> rho_th count_th NN_check <triangle_path> <mf_path>\n", argv[0]);
    return EXIT_FAILURE;
  }

//argv[0]=./main; argv[1]=<inputfilename>; argv[2]=ro_th; argv[3]=LL; argv[4]=<triangle_path>; argv[5]= <outfile MF> argc=6

	//printf(" \n %s \t %s \t %s  \n",argv[0],argv[1],argv[2]);
//return EXIT_SUCCESS;

  if(access(argv[1], F_OK)!=0)     {
    printf("File %s does not exist\n Now exiting to system ....\n", argv[1]);
    return EXIT_FAILURE; } //checks the <inputfilename> is there or not
  
 // ReadFile(argv[1], argv[2], argv[3], argv[4]);
  iso= atof(argv[2]); //INDEX=iso; //INDEX=rho_th; atof reads the string and convert to float
  count_th=atoi(argv[3]); LL=1.0;
  NN_check=atoi(argv[4]);
  strcpy(fout, argv[5]);
  strcpy(fname, argv[6]);
  printf("rho_th=%lf, count_th=%d \n.......................................................\n",iso, count_th);
  //printsize();
  //Read_density(); //read the density file
  //printsize();
  printf("Read_density\n");
  
  //NN_check=-1;
  
  NC = Drive_cluster(argv[1]);
  //printsize();
  //printf("NC = %d\n", NC);  
  /*storing Minkowski Functionals into file */
  
  /*outfile = fopen(argv[5], "a");
  for(ii=0; ii<NC; ii++)    {
    PrepareGrid(ii+1);
    run();
    printsize();
    index=store_triangle(NC); 

    area=compute_area();

    vol=compute_vol();
    //imc=compute_imc();

    imc = compute_curvature();

    genus=compute_genus();

    if(vol>0)
      Shape1 = 3.*(vol / area); //Shape1 =T
    else
      Shape1 = -3.*(vol / area);
    Shape2 = area / imc;  //Shape2=B
    if(genus)
      Shape3 = imc/(4.0 * M_PI * abs(genus)); //Shape3=L
    else
      Shape3 = imc/(4.0 * M_PI);

    if(Shape1 > Shape2)
      {
        s1 = Shape2;
        s2 = Shape1;
      }
    else
      {
        s1 = Shape1;
        s2 = Shape2;
      }
    if(s2 > Shape3){
      s3 = s2;
      s2 = Shape3;
    }
    else
      s3 = Shape3;
    
    if(s1 > s2)
      {
        temp = s2;
        s2 = s1;
        s1 = temp;
      }
    
    Planarity = (s2 - s1) / (s2 + s1);
    Filamentarity = (s3 - s2) / (s3 + s2);
    
    if(index>0)
      fprintf(outfile,"INDEX= %f\t index= %d\t area= %lf\t vol= %lf\t imc= %lf"
              "\t genus= %lf\t Shape1= %lf\t Shape2= %lf\t Shape3= %lf\t Planarity= %lf\t Filamentarity= %lf\n", INDEX, index, area, vol, imc,
              genus, Shape1, Shape2, Shape3, Planarity, Filamentarity);
    clean_all();
    printf("\n INDEX= %f\t index= %d\t area= %lf\t vol= %lf\t imc= %lf"
              "\t genus= %lf\t Shape1= %lf\t Shape2= %lf\t Shape3= %lf\t Planarity= %lf\t Filamentarity= %lf\n \n", INDEX, index, area, vol, imc,
              genus, Shape1, Shape2, Shape3, Planarity, Filamentarity);
  }
  printf("NC = %d\n", NC); 
  fclose(outfile);*/
  
//satadru();
  return EXIT_SUCCESS;
}
