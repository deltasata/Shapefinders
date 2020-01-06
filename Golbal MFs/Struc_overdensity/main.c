# include "MarchingCube.h"
# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <float.h>
# include <stdbool.h>
#include <mpi.h>

int NI, NJ, NK, NIJK, NIJK1;
double * rho;
char fin[300];

FILE *fp2;



int ncpoints, count1; //ncpoints is sum of all points (rho>=rho_th) counted in DFSinteractive(). count1 is counted as rho>=rho_th_min



/*# define get_x(index) (index/(NJ*NK))
# define get_y(index) ((index/NK)%NK)
# define get_z(index) (index%NK)*/

int get_x(int index){return index/(NK*NJ);}
int get_y(int index){return (index/NK)%NJ;}
int get_z(int index){return index%NK;}

void get_ijk(int index, int ijk[3]){ijk[0]=index/(NK*NJ); ijk[1]=(index/NK)%NJ;  ijk[2]=index%NK;}

int get_index(int i, int j, int k){return ((i*NJ + j)*NK + k);}

# define clamp(a, amin, amax) ((a)>(amax)?(amax):((a) < (amin)? (amin):(a)))
void min_max(int a, int *min, int *max){if(a>*max){*max=a;} if(a<*min){*min=a;}}

void DataFeed_SFcal(double, int);

//========================================= main starts ==============================================================

int main(int argc, char **argv)  {
  int unused __attribute__((unused));
  //FILE *outfile;
  //int NT=0;
  
  if(argc!=7) {
    printf("usage:\n\t %s <inputfilename> rho_th_min rho_th_max step <triangle_path> <mf_path>\n", argv[0]);
    return EXIT_FAILURE;
  }


  if(access(argv[1], F_OK)!=0)     {
    printf("File %s does not exist\n Now exiting to system ....\n", argv[1]);
    return EXIT_FAILURE; } //checks the <inputfilename> is there or not
  
  double rho_th_min, rho_th_max; int step;
 // ReadFile(argv[1], argv[2], argv[3], argv[4]);
 //INDEX=iso; //INDEX=rho_th; atof reads the string and converts to float
  rho_th_min= atof(argv[2]); rho_th_max= atof(argv[3]); step=atoi(argv[4]);
  LL=1.0; count_th=1;
  strcpy(fout, argv[5]); //for saving the triangles
  strcpy(fname, argv[6]);


	int ierr,my_id,an_id,num_procs,root_process=0,num_rows_to_receive; int tag=1;
        MPI_Status status;
        ierr = MPI_Init(&argc, &argv);
        
        ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
        ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	int last_id=num_procs-1;
  
  NN_check=-1;
  
  strcpy(fin,argv[1]);
  
  FILE *fp;
  if ((fp = fopen(fin, "r")) == NULL)    
  {
    fprintf(stderr,"FILE Cannot open %s...proccess no.=%d \n", fin, my_id);
    exit(1);
  }
  
  int index;
  int ii, jj, kk,  NC=1; //NC does nothing
  double ff, rhob, dmin, dmax, sig;
   
  double temp1,temp2; //of no use
  fscanf(fp, "#%d%d%d%lf%lf", &NI, &NJ, &NK,&temp1,&temp2); //just to match the input file we use temp1 and temp2
  _size_x=NI+1; _size_y=NJ+1; _size_z=NK+1;
  NIJK = (int)(NI*NJ*NK); NIJK1=(NI+1)*(NJ+1)*(NK+1); // needed to store the final data[] array
  
  rho = (double*) calloc(NIJK, sizeof(double));
  data = (double*) calloc(_size_x*_size_y*_size_z, sizeof(double));
  
  

  rhob=0; sig=0; dmin=FLT_MAX; dmax=-1.*FLT_MAX;

  count1=0;

  for(ii=0; ii<NI; ii++)
    for(jj=0; jj<NJ; jj++)
      for(kk=0; kk<NK; kk++)   {

	index = (ii*NJ + jj)*NK + kk;
	fscanf(fp, "%lf", &rho[index]); if(rho[index]>=rho_th_min){count1++;}
	rhob += rho[index];
        sig += (rho[index]*rho[index]);
        dmin = (dmin < rho[index]) ? dmin : rho[index];
        dmax = (dmax > rho[index]) ? dmax : rho[index];
      }
  fclose(fp);
  
  for(ii=0;ii<NIJK1;ii++){data[ii]=-1.0;}


  rhob /=(NI*NJ*NK);
  sig = sqrt((sig/(NI*NJ*NK)) - (rhob*rhob));
  
   sprintf(fshape,"%sShapefinders_iso%lf_to_%lf",fname, rho_th_min,rho_th_max);
  //for deleting data from existing file, if any
	if(my_id==last_id)
	{
		printf("....................... \n Read_density done............\n");
		printf(" The input data box: NI= %d, NJ= %d, NK= %d, NIJK= %d \n",NI,NJ,NK,NIJK);
  		printf("rho_th_min=%lf, rho_th_max=%lf, step=%d, count1 for rho_th_min=%d \n", rho_th_min, rho_th_max, step, count1);
  		printf("*********************************************************************\n\n");
  		
  		fp2=fopen(fshape, "w");
  		fclose(fp2);
				
	}

  //double rho_increment=1.0; //linear spacing
  double rho_increment=pow((rho_th_max/rho_th_min),(1.0/(step-1.0))); //logspace
  
  //double start=rho_th_min+my_id*rho_increment; double last_rh_th=rho_th_min+(last_id-1.0)*rho_increment; //linear spacing

  double start=rho_th_min*pow(rho_increment,(my_id)); double last_rh_th=rho_th_min*pow(rho_increment,(last_id-1)); //logspace

  double rho_th;
  int working_procs=num_procs-1;
  if(my_id == last_id) //last_id is the clerical processor
  {
        int i_working_procs=0; 
  	bool all_done=false;
  	while(all_done==false)
	{
	        an_id=-1; //not using it
		ierr = MPI_Recv( &an_id, 1, MPI_INT, MPI_ANY_SOURCE,tag, MPI_COMM_WORLD, &status);
		last_rh_th=last_rh_th*rho_increment; //logspace
		//last_rh_th=last_rh_th+rho_increment; //for linespacing
		
		if(last_rh_th<rho_increment*rho_th_max)
		{
			if(last_rh_th>rho_th_max){last_rh_th=rho_th_max;}
		        ierr = MPI_Send(&last_rh_th, 1, MPI_DOUBLE, an_id, tag, MPI_COMM_WORLD);
		}
		else
		{
		        double dummy_rh=-1.0;
		        ierr = MPI_Send( &dummy_rh, 1, MPI_DOUBLE, an_id, tag, MPI_COMM_WORLD);
		        i_working_procs++;
		        if(i_working_procs==working_procs)
		        {all_done=true;}
		
		}
	
	}
	printf("Last proc end:\n");
	
  
  }
  else
  {
        rho_th=start;
  
  	while(rho_th>=0)
  	{
  		
  		DataFeed_SFcal(rho_th, my_id);
		
                ierr = MPI_Send( &my_id, 1, MPI_INT, last_id, tag, MPI_COMM_WORLD);
  	        ierr = MPI_Recv( &rho_th, 1, MPI_DOUBLE, MPI_ANY_SOURCE,tag, MPI_COMM_WORLD, &status);
  	
  	}
  
  	
  
  }
  
  ierr = MPI_Finalize();	
  return EXIT_SUCCESS;
}


//=================================================================================================================================
//========================================= main ends =============================================================================

void DataFeed_SFcal(double rho_th, int my_id)
{
	int NN=1; bad_NE=0; //NN does nothing
	//printf("....starting rho_th=%.16lf....bad_NE=%d \n", rho_th, bad_NE);
	iso=rho_th;
	ncpoints=0;
	int i,j,k,in,jn,kn, index_old;
  	int indi; indi=0;
  	for(in=0;in<_size_x;in++)
  	{
  		for(jn=0;jn<_size_y;jn++)
  		{
  			for(kn=0;kn<_size_z;kn++)
  			{
  				i=in%NI; j=jn%NJ; k=kn%NK;
  				index_old=get_index(i,j,k); 

  				data[indi]=rho[index_old];
  				if((i==0 || j==0 || k==0))
  				{ 
  					if(data[indi]>=rho_th) {data[indi]=0.9*rho_th;}
  				}
				if(data[indi]>=rho_th){ncpoints++;}
				indi++;
  			}
  			
  			
  		}
  	}

   
//********************************* running Merchingcube ****************************************************************************
//===================================================================================================================================
		double imc_used;
		
		//init_all();
	
		run(NN); //now everything this calculated in run() and outputs are saved into global variables declared in Marchingcube.h
		
		
		/*no_of_tri=store_triangle(NN); 

    		area=compute_area();

    		vol=compute_vol();
    		//imc=compute_imc();

    		imc = compute_curvature();

    		genus=compute_genus_imc();*/   //genus=compute_genus() for cal_genus.c
    		
    		if(imc2!=imc2){imc_used=imc;} else{imc_used=imc2;}
    		
    		double Shape1=3.0*fabs(vol)/area;
    		double Shape2=area/imc_used;
    		double Shape3=imc_used/(4.0*M_PI/**(1.0+fabs(genus))*/);
    		double Shape31=imc_used/(4.0*M_PI*(1.0+fabs(genus)));
    		double T=Shape1; double B=Shape2; double L=Shape3; double L1=Shape31;
    		
    		char *text;
    		if(L>=B && B>=T){text="L>=B>=T";} 
    		else { text="L,B,T not in order"; }
    		double P=(B-T)/(B+T); double F=(L-B)/(L+B); double F1=(L1-B)/(L1+B);
    
    		printf("\n********************************* rho_th=%lf, my_id=%d ***************************************************\n <<<<<<<<<<<<<<<<<<<< diff_NE=%lf, bad_NE=%d, zero_tri=%d >>>>>>>>>>>>>>>>>>>> \niso=%.16lf, ncpoints=%d,   Vol=%.16lf,   area=%lf, imc=%lf, imc1=%lf, imc2=%lf, genus=%lf, T=%lf, B=%lf, L=%lf, Planarity=%lf, Filamentarity=%lf, Filamentarity1=%lf \n %s \n==================================================================== \n",rho_th,my_id,diff_NE, bad_NE,zero_tri,iso,ncpoints,vol,area,imc,imc1,imc2,genus,T,B,L,P,F,F1,text);
  		//printf("\n \n %lf \t  %lf \t %lf \n \n ",vol,area,imc2);
  		fp2 = fopen(fshape, "a");
  		fprintf(fp2,"%.16lf \t %d \t %.16lf \t %.16lf \t %.16lf \t %.16lf \t %.16lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \n", iso,ncpoints,vol,area,imc,imc1,imc2,genus,T,B,L,P,F,F1);
  		fclose(fp2);



}

