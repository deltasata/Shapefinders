//this code works correctly when the periodicity is excluded in the data box fed to it. 0 \equiv NI (the box will always have data from 0 to NI-1). If the box has periodic data (i.e. 0 \eqiv NI-1), to make this correct, one needs to change iPeriod and the final arrangement part: %NI should be replaced to NI-1. Also in the final part we do not need a box six=ze of NI+1, instead till NI. 
//Author=Satadru Bag

//compile: mpicc -o r1.out NC_finder_all_parallel_new_struc.c -lm
//run: mpirun -np 64 ./r1.out
# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <float.h>
# include <stdbool.h>
#include <mpi.h>
//int N;

int NI, NJ, NK, NIJK, NIJK1;
double * rho, rho_th;

//================ input and output files ==============================================================
char fin[300]="/home/aai/satadru/Juhan_data/NEW_1mpc/density_scl_m_230419_dr12_a12_d1.txt"; //input (density) field file
char fout[300]="struc.dat"; // output file
//output file structure: rho_th,NC,ff, percolation_no, NN_max, count_max,ncpoints,count1, LCS

//======================================================================================================

double rho_min=10000.0; double rho_max=0.0;
double rho_min1=10000.0;

int NN_max;
int count_max; 


int par[3];
int percolation_no;

int _size_x,_size_y,_size_z;



int ncpoints, count1,not; //ncpoints is sum of all points counted in DFSinteractive(). count1 is counted as rho>=rho_th

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

int ijk_shift[3]={0,0,0}; //shift of origin to arrange each cluster (fragmented or not fragmented)

typedef struct Stack {
  int number;  //stores the index of a point
  int parity[3]; //stores xp, yp, zp: need to identify the shift of (0,0,0) while rearranging a fragmented cluster
  struct Stack *next;
}*Stack;


void push(struct Stack**, int*, int);
int pop(struct Stack**,int*);
int DFSinteractive(int*, int);
int iPeriod(int*,int *, int);
int DriveClusterFinder( int *, double,  double*);
void mapping(int*, double*, double, int, double*);




int main ( int argc, char *argv[] )
{
  if(argc!=4)
  {
	printf("usage:\n\t %s <rho_th_min> <rho_th_max> <no. of steps>\n", argv[0]);
	return EXIT_FAILURE;
  }
  double rh_min=atof(argv[1]); double rh_max=atof(argv[2]); int Step=atoi(argv[3]);

  
  FILE *fp; FILE *fp2;
  if ((fp = fopen(fin, "r")) == NULL)    
  {
    fprintf(stderr,"FILE Cannot open %s \n", fin);
    exit(1);
  }
  
  int index;
  int ii, jj, kk, *cluster,  NC;
  //int *xmin, *xmax, *ymin, *ymax, *zmin, *zmax;
  double ff, ff1, rhob, dmin, dmax, sig, LCS;
 
  int unused __attribute__((unused));


	int ierr,my_id,an_id,num_procs,root_process=0,num_rows_to_receive; int tag=1;
        MPI_Status status;
        ierr = MPI_Init(&argc, &argv);
        
        ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
        ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	int last_id=num_procs-1;

//for deleting data from existing file, if any
	if(my_id==last_id)
	{
		fp2 = fopen(fout, "w");
  		fclose(fp2);
	}
  
  //printf("rho_th=%lf \n", rho_th);
  
  double temp1,temp2; //of no use
  fscanf(fp, "#%d%d%d%lf%lf", &NI, &NJ, &NK,&temp1,&temp2); //just to match the input file we use temp1 and temp2
  _size_x=NI+1; _size_y=NJ+1; _size_z=NK+1;
  NIJK = (int)(NI*NJ*NK); NIJK1=(NI+1)*(NJ+1)*(NK+1); // needed to store the final data[] array
  

  
  
  //for linear spacing
  double rho_increment=(rh_max-rh_min)/(Step-1.0);
  double start=rh_min+my_id*rho_increment; 
  double last_rh_th=rh_min+(last_id-1.0)*rho_increment;


  //for log spacing
  //double rho_increment=pow((rh_max/rh_min),(1.0/(Step-1.0))); //logspace
  //double start=rh_min*pow(rho_increment,(my_id)); 
  //double last_rh_th=rh_min*pow(rho_increment,(last_id-1));
  
  double rho_th;
  int working_procs=num_procs-1;
  if(my_id == last_id) //last_id is the clerical processor
  {
    	int i_working_procs=0; 
  	bool all_done=false;
  	while(all_done==false)
	{
	        an_id=-1;
		ierr = MPI_Recv( &an_id, 1, MPI_INT, MPI_ANY_SOURCE,tag, MPI_COMM_WORLD, &status);
		//last_rh_th=last_rh_th*rho_increment; //logspace
		last_rh_th=last_rh_th+rho_increment; //for linear spacing
		
		if(last_rh_th<=rh_max)
		{
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
  
	rho = (double*) calloc(NIJK, sizeof(double));
	rhob=0; sig=0; dmin=FLT_MAX; dmax=-1.*FLT_MAX;

	for(ii=0; ii<NI; ii++)
		for(jj=0; jj<NJ; jj++)
			for(kk=0; kk<NK; kk++)
			{
				index = (ii*NJ + jj)*NK + kk;
				fscanf(fp, "%lf", &rho[index]);
				if(rho[index]>rho_max){rho_max=rho[index];}
				if(rho[index]<rho_min){rho_min=rho[index];}
				if(rho[index]<rho_min1 && rho[index]!=0){rho_min1=rho[index];}
				rhob += rho[index];
				sig += (rho[index]*rho[index]);
				dmin = (dmin < rho[index]) ? dmin : rho[index];
				dmax = (dmax > rho[index]) ? dmax : rho[index];
			}
  	fclose(fp);
	//rhob /=(NI*NJ*NK);
	//sig = sqrt((sig/(NI*NJ*NK)) - (rhob*rhob));
	cluster = (int *) calloc(NIJK, sizeof(int));
        
        rho_th=start;
  
  	while(rho_th>=0)
  	{
  		count_max=0; ncpoints=0; count1=0; percolation_no=0;	
  		NC = DriveClusterFinder(cluster, rho_th, &ff);
  		ff1=((double) ncpoints)/((double) NIJK);
  		LCS=((double) count_max)/((double) ncpoints);
  		//printf("rho_th=%lf, NC = %d; ff=%lf; percolation_no=%d; NN_max=%d, count_max=%d  \n", rho_th, 	NC,ff,percolation_no,NN_max,count_max);

		printf("%.16lf \t %d \t %lf \t %d \t %d \t %d \t %d \t %d \t %lf \n", rho_th,NC,ff, percolation_no, NN_max, count_max,ncpoints,count1,LCS);
  	
		fp2 = fopen(fout, "a");
  		fprintf(fp2, "%.16lf \t %d \t %lf \t %d \t %d \t %d \t %d \t %d \t %.16lf \n", rho_th,NC,ff, percolation_no, NN_max, count_max,ncpoints,count1, LCS);
  		fclose(fp2);

                ierr = MPI_Send( &my_id, 1, MPI_INT, last_id, tag, MPI_COMM_WORLD);
  	        ierr = MPI_Recv( &rho_th, 1, MPI_DOUBLE, MPI_ANY_SOURCE,tag, MPI_COMM_WORLD, &status);
  	
  	}
  
  	
  
  }
  

  //printf("NC = %d; percolation_no=%d; NN_max=%d, count_max=%d  \n", NC,percolation_no,NN_max,count_max);
    

 
   ierr = MPI_Finalize();
  
   
  
  //printf("count1=%d, ncpoints=%d : they should match. NIJK=%d: ff=%lf, ff1=%lf \n",count1, ncpoints,NIJK,ff,ff1);
   
}

//============================ main ends ===================================================================================================
/*********************************** push a element v in the stack S **********************************************************************/
void push(struct Stack **start, int v[5], int x) //x=1: i direction. x=2: j direction . x=3: k direction.
{
	if(v[4]>-1){par[x-1]=-1;ijk_shift[x-1]=0;}

	int ijk[3]; get_ijk(v[0],ijk);
  if(v[x]==par[x-1] && ijk[x-1]>ijk_shift[x-1]) {ijk_shift[x-1]=ijk[x-1];} //only to check xth parity. If it is 1, then only check if xth cordinate is maximum
    
  struct Stack *new;
  //printf("v = %ld\n", v);
  new = (struct Stack*)malloc(sizeof(struct Stack));
  new->number=v[0];
  new->parity[0]=v[1]; new->parity[1]=v[2]; new->parity[2]=v[3]; // store the index and parity of that point
  new->next=*start;
  *start=new;
}

/* pop an element from the stack */ //returns the index and gives parities of the point in v[] 
int pop(struct Stack **start, int v[4]) {
  struct Stack *temp;
  int value;
  if(*start==NULL) 
    return -1;
  else {
    value=v[0]=(*start)->number; v[1]=(*start)->parity[0]; v[2]=(*start)->parity[1]; v[3]=(*start)->parity[2];
    temp=(*start);
    (*start)=(*start)->next;
    free(temp);
  }
  return value;
}
//********************************************************************************************************************************
//========== for periodic boundary condition =====================================================================================

int iPeriod(int nb[4],int nb_dummy[5], int n)
{
	int p;
	int I=nb[0]; int k; for(k=0;k<4;k++){nb_dummy[k]=nb[k];} nb_dummy[4]=-1;
	int ijk[3]; int NIJK_array[3]={NI,NJ,NK}; //ijk[3] stores i,j,k for a given index I
	//ijk[0]=get_x(I); ijk[1]=get_y(I); ijk[2]=get_z(I);
	get_ijk(I,ijk);
	//p=1: parity of the points with the point found first. p=-1 connected points but separated by the wall. 
	if((ijk[abs(n)-1]+abs(n)/n)<0 || (ijk[abs(n)-1]+abs(n)/n)>=NIJK_array[abs(n)-1]) 
	{
		p=-1;
		//to prevent a special case: 
		if((ijk[abs(n)-1]+abs(n)/n)>=NIJK_array[abs(n)-1] && nb[abs(n)]==1 && par[abs(n)-1]==1){nb_dummy[4]=abs(n)-1;}
	}
	else{p=1; }
	
	nb_dummy[abs(n)]=p*nb_dummy[abs(n)]; // stores x,y,z parity of the point. Note only the ith (x or y or z) coordinate parity may change when proceed in ith direction
	
	
	ijk[abs(n)-1]=(ijk[abs(n)-1]+abs(n)/n+NIJK_array[abs(n)-1])%NIJK_array[abs(n)-1]; //only the required element of i,j,k is changed with pbc
	//printf(" I=%d, 1st=%d, n=%d, I+n=%d, 2nd= %d \n",I, get_x(I),n,I+n, ijk[0]);
	
	if(ijk[0]<0 || ijk[1]<0 || ijk[2]<0){printf("iPeriod PROBLEM: \n \n");} //check this..not sure
	 
	  I=(ijk[0]*NJ+ijk[1])*NK+ijk[2]; 
	  nb_dummy[0]=I; // stores the correct index I of the point to nb_dummy[0]; not to nb[0] because nb[] represent the mother point, not this point.
	  return I;
	
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//>>>>>>>>>>>>>>>>>>>>>>>>>>> DriveClusterFinder starts <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
int DriveClusterFinder(int *cluster, double rho_th, double *ff) {
  int tag = -1, flag=1, NC=0;
 // int ii, jj, kk;
  mapping(cluster, rho, rho_th, tag, ff);
  

//enters with flag=1
  while(flag>=0) {
    //printf("flag = %d\n", flag);
    flag = DFSinteractive(cluster, flag); //returns flag=-1 if not any other point found, while loop breaks
    if(flag>0) NC = flag-1;
  }
  
  return NC;			/* no of clusters */
}


void mapping(int *cluster, double *rho, double rho_th, int tag, double *ff) {

  int index;
  int ii, jj, kk;
  //NIJK = NI*NJ*NK;
  
  for (ii=0; ii<NI; ii++) 
    for (jj=0; jj<NJ; jj++) 
      for (kk=0; kk<NK; kk++) {
	index = ii*NJ*NK + jj*NK + kk;

	if(rho[index]>=rho_th) {
	  cluster[index] = tag;
	  count1++;
	}
	else
	  cluster[index] = 0;
      }
  //printf("total count = %ld \t %ld\n", count1, NIJK);
  *ff = (double)count1/(double)NIJK;
}

/*int check_percolation(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax)
{
	
	printf("%d, %d, %d ::%d, %d, %d, %d, %d, %d \n",NI,NJ,NK,xmin,xmax,ymin,ymax,zmin,zmax);
	if(abs(xmax-xmin)>=(NI-1) || abs(ymax-ymin)>=(NJ-1) || abs(zmax-zmin)>=(NK-1))
	{return 1;}
	else {return -1;}	
}*/


//=========================================================================================================================================
int DFSinteractive(int *cluster, int NN) 
{

	/*if(NN>N-1)
	{
		printf("...NN>N-1.. N=%d, NN=%d \n",N,NN);
		N=NN*2;
		xmin=(int *) realloc(xmin, N*sizeof(int));xmax=(int *) realloc(xmax, N*sizeof(int));ymin=(int *) realloc(ymin, N*sizeof(int));ymax=(int *) realloc(ymax, N*sizeof(int));zmin=(int *) realloc(zmin, N*sizeof(int));zmax=(int *) realloc(zmax, N*sizeof(int));
	
	}*/

	int count=0;
	

  struct Stack *S=NULL;
  int ii=0, vv;  
  //int NIJK = NI*NJ*NK;
  
  while((cluster[ii]>-1)&&(ii<NIJK)) {
    //printf("ok \n");
    ii++;
  }
  /* return -1 if there are no clusters in the distributions */
  //printf("%ld\n", ii);
  //printf("%d\n", cluster[ii]);

	//printf("Satadru: %d \n",NN);

//if(NN==1){printf("start check cluster: %d , %d, %d, %d \n",cluster[(0*NJ+40)*NK+40],cluster[(100*NJ+200)*NK+100],cluster[(100*NJ+100)*NK+0],cluster[(100*NJ+100)*NK+200]);}

  if(ii>=NI*NJ*NK)
    return -1;
    
  cluster[ii]=NN; count++; //assummes parity[3]=1 for all indices
  //zmax[NN] =zmin[NN]=ijk_shift[2]= get_z(ii); ymax[NN] =ymin[NN]= ijk_shift[1]=get_y(ii);xmax[NN] =xmin[NN]= ijk_shift[0]=get_x(ii);
  
  ijk_shift[2]= get_z(ii);   ijk_shift[1]=get_y(ii);  ijk_shift[0]=get_x(ii); 
  
  //printf("ii = %ld\n", ii);

//the convention is: 1= move 1 in +x dir, -1= move 1 in -x dir
//                   2= move 1 in +y dir, -2= move 1 in -y dir
//                   3= move 1 in +z dir, -3= move 1 in -z dir  

  int nb[4]={ii,1,1,1};par[0]=1;par[1]=1;par[2]=1;//nb[0] stores index. nb[1,2,3]=xp,yp,zp store parity
  int nb_dummy[5]={ii,1,1,1,-1}; int nxt;

//iPeriod changes nb_dummy[] according to the index and parity. Also sets nb_dummy[0]=nxt.  nb[] is not changed, stores info of the mother point (index ii)
  nxt = iPeriod (nb,nb_dummy, 3); //checks the next point in +z direction
  if((nxt!=ii)&&(cluster[nxt] == -1)) {push(&S, nb_dummy,3);cluster[nxt]=NN;count++;} 
  

  nxt = iPeriod (nb,nb_dummy, -3);
  if((nxt!=ii)&&(cluster[nxt] == -1)) {push(&S, nb_dummy,3);cluster[nxt]=NN;count++;}
 

  nxt = iPeriod (nb,nb_dummy, 2);//checks the next point in +y direction
  if((nxt!=ii)&&(cluster[nxt] == -1)) {push(&S, nb_dummy,2);cluster[nxt]=NN;count++;}
  

  nxt = iPeriod (nb,nb_dummy, -2);
  if((nxt!=ii)&&(cluster[nxt] == -1)) {push(&S, nb_dummy,2);cluster[nxt]=NN; count++;}
  

  nxt = iPeriod (nb,nb_dummy, 1);//checks the next point in +x direction
  if((nxt!=ii)&&(cluster[nxt] == -1)) {push(&S, nb_dummy,1);cluster[nxt]=NN;count++;}
  
  
  nxt = iPeriod (nb,nb_dummy, -1);
  if((nxt!=ii)&&(cluster[nxt] == -1)) {push(&S, nb_dummy,1);cluster[nxt]=NN;count++; }
  
    

    
/********** repeat the procedure for each point stored in the stack until it is empty **************************************************/

  while(S!=NULL) 
  {
  	
    	//printf("Entered the last while loop \n");
    	vv = pop(&S,nb);
    	if (vv==-1) {return NN+1;printf("CHeck: weired....\n \n \n ");} //did not understand properly, does not enter here

  	nxt = iPeriod (nb,nb_dummy, 3); 
  	if((nxt!=vv)&&(cluster[nxt] == -1)) {push(&S, nb_dummy,3);cluster[nxt]=NN;count++;} 
  

  	nxt = iPeriod (nb,nb_dummy, -3);
  	if((nxt!=vv)&&(cluster[nxt] == -1)) {push(&S, nb_dummy,3);cluster[nxt]=NN;count++;}
 

  	nxt = iPeriod (nb,nb_dummy, 2);
  	if((nxt!=vv)&&(cluster[nxt] == -1)) {push(&S, nb_dummy,2);cluster[nxt]=NN;count++;}
  

  	nxt = iPeriod (nb,nb_dummy, -2);
  	if((nxt!=vv)&&(cluster[nxt] == -1)) {push(&S, nb_dummy,2);cluster[nxt]=NN; count++;}
  

  	nxt = iPeriod (nb,nb_dummy, 1);
  	if((nxt!=vv)&&(cluster[nxt] == -1)) {push(&S, nb_dummy,1);cluster[nxt]=NN;count++;}
  
  
  	nxt = iPeriod (nb,nb_dummy, -1);
  	if((nxt!=vv)&&(cluster[nxt] == -1)) {push(&S, nb_dummy,1);cluster[nxt]=NN;count++;}
      

  }
  if(count>count_max){count_max=count; NN_max=NN;}
 ncpoints=ncpoints+count;	

 	
 	//printf("NC=%d: count= %d :",NN,count); 
  /*if(check_percolation(xmin[NN], xmax[NN], ymin[NN], ymax[NN], zmin[NN], zmax[NN])>0)
  {
  	printf("Percolation happend. NC= %d \n", NN);
  }*/
  
  //printf("NC= %d:: ijk_shift= {%d, %d, %d}...stroing points... \n",NN,ijk_shift[0],ijk_shift[1],ijk_shift[2]);
  
  ijk_shift[0]=(ijk_shift[0]+1)%NI; ijk_shift[1]=(ijk_shift[1]+1)%NJ; ijk_shift[2]=(ijk_shift[2]+1)%NK;
  
   
  
  int in,jn,kn,i,j,k, index_old, index_new; bool clause =true;
  
  //int indi;
  //indi=0;
  
  
  
 if(count>=400)
 {
  	/*for(in=0;in<_size_x;in++)
 	{
  		for(jn=0;jn<_size_y;jn++)
  		{
  			for(kn=0;kn<_size_z;kn++)
  			{
  				i=(in+ijk_shift[0])%NI; j=(jn+ijk_shift[1])%NJ; k=(kn+ijk_shift[2])%NK;
  				index_old=get_index(i,j,k); index_new=get_index(in,jn,kn);
  			
  			
  				if(cluster[index_old]==NN)
  				{
  					if((in==0 || jn==0 || kn==0))
  					{
  						percolation_no++;
  						return NN+1;
  					}	
  			
  				}
  			
  			}
  		}
    	}*/

	in=0;
  	for(jn=0;jn<_size_y;jn++)
  	{
  		for(kn=0;kn<_size_z;kn++)
  		{
  			i=(in+ijk_shift[0])%NI; j=(jn+ijk_shift[1])%NJ; k=(kn+ijk_shift[2])%NK;
  			index_old=get_index(i,j,k); index_new=get_index(in,jn,kn);
  			
  			
  			if(cluster[index_old]==NN)
  			{
  				percolation_no++;
  				return NN+1;	
  			
  			}
  			
  		}
  	}
  	
  	jn=0;
  	for(in=0;in<_size_x;in++)
  	{
  		for(kn=0;kn<_size_z;kn++)
  		{
  			i=(in+ijk_shift[0])%NI; j=(jn+ijk_shift[1])%NJ; k=(kn+ijk_shift[2])%NK;
  			index_old=get_index(i,j,k); index_new=get_index(in,jn,kn);
  			
  			
  			if(cluster[index_old]==NN)
  			{
  				percolation_no++;
  				return NN+1;	
  			
  			}
  			
  		}
  	}
  	
  	kn=0;
  	for(in=0;in<_size_x;in++)
  	{
  		for(jn=0;jn<_size_y;jn++)
  		{
  			i=(in+ijk_shift[0])%NI; j=(jn+ijk_shift[1])%NJ; k=(kn+ijk_shift[2])%NK;
  			index_old=get_index(i,j,k); index_new=get_index(in,jn,kn);
  			
  			
  			if(cluster[index_old]==NN)
  			{
  				percolation_no++;
  				return NN+1;	
  			
  			}
  			
  		}
  	}
  
  }		
  return NN+1;
}
