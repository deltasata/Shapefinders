//this code works correctly when the periodicity is excluded in the data box fed to it. 0 \equiv NI (the box will always have data from 0 to NI-1). If the box has periodic data (i.e. 0 \eqiv NI-1), to make this correct, one needs to change iPeriod and the final arrangement part: %NI should be replaced to NI-1. Also in the final part we do not need a box size of NI+1, instead till NI. 

# include "MarchingCube.h"
# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <float.h>
# include <stdbool.h>
//int N;

int NI, NJ, NK, NIJK, NIJK1;
double * rho, rho_th;
char fin[300];
int par[3];

FILE *fp2;
FILE *fp3;



int ncpoints, count1; //ncpoints is sum of all points counted in DFSinteractive(). count1 is counted as rho>=rho_th
double tot_vol;

int percolation_no;
int count_max, NN_max; //which cluster has maximum points inside it.
double count_max_vol; //volume of that cluster which has max count; can not be calculated always, if count< threshold, not goes to triangulation

double vol_max;
int NN_max_vol;
int vol_max_count; //count of that cluster which has max volume

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


double rho_0=0.001; //for the current Durham problem

double void_func(double x){return rho_0/(0.000001+x);} //suitable when x (rth) is of the order unity?

int Drive_cluster(char *f1)
{

	strcpy(fin, f1);
  
  FILE *fp;
  if ((fp = fopen(fin, "r")) == NULL)    
  {
    fprintf(stderr,"FILE Cannot open %s \n", fin);
    exit(1);
  }
  
  int index;
  int ii, jj, kk, *cluster,  NC;
  //int *xmin, *xmax, *ymin, *ymax, *zmin, *zmax;
  double ff, rhob, dmin, dmax, sig;
 
  int unused __attribute__((unused));

  double rho_th1=iso;
  rth=rho_th1; //a global variable set here (can be used in other files)

  rho_th = void_func(iso);
  iso=void_func(iso);

  
  double temp1,temp2; //of no use
  fscanf(fp, "#%d%d%d%lf%lf", &NI, &NJ, &NK,&temp1,&temp2); //just to match the input file we use temp1 and temp2
  _size_x=NI+1; _size_y=NJ+1; _size_z=NK+1; //why _size_x is NI+1, not just NI?
  NIJK = (int)(NI*NJ*NK); NIJK1=(NI+1)*(NJ+1)*(NK+1); // needed to store the final data[] array
  
  rho = (double*) calloc(NIJK, sizeof(double));
  data = (double*) calloc(_size_x*_size_y*_size_z, sizeof(double));
  
  

  rhob=0; sig=0; dmin=FLT_MAX; dmax=-1.*FLT_MAX;

  for(ii=0; ii<NI; ii++)
    for(jj=0; jj<NJ; jj++)
      for(kk=0; kk<NK; kk++)   {

	index = (ii*NJ + jj)*NK + kk;
	fscanf(fp, "%lf", &rho[index]);

	rho[index]=void_func(rho[index]);
	rhob += rho[index];
        sig += (rho[index]*rho[index]);
        dmin = (dmin < rho[index]) ? dmin : rho[index];
        dmax = (dmax > rho[index]) ? dmax : rho[index];
      }
  fclose(fp);
  
  for(ii=0;ii<NIJK1;ii++){data[ii]=-1.0;}


  rhob /=(NI*NJ*NK);
  sig = sqrt((sig/(NI*NJ*NK)) - (rhob*rhob));
  
 // N=20000; 
  ncpoints=0; count1=0; tot_vol=0;
  percolation_no=0; count_max=0;  NN_max=0; NN_max_vol=0; vol_max=0.0; count_max_vol=-1; vol_max_count=0; bad_NE=0;
 // xmin=(int *) calloc(N, sizeof(int));xmax=(int *) calloc(N, sizeof(int));ymin=(int *) calloc(N, sizeof(int));ymax=(int *) calloc(N, sizeof(int));zmin=(int *) calloc(N, sizeof(int));zmax=(int *) calloc(N, sizeof(int));
  
  printf(" The input data box: NI= %d, NJ= %d, NK= %d, NIJK= %d \n",NI,NJ,NK,NIJK);

  printf("**********************************************\n\n");

  /* printf("\t mean : %e\t sigma : %e\n", rhob, sig); */
  /* printf("\t dmin : %e\t  dmax : %e\n", dmin, dmax); */

  cluster = (int *) calloc(NIJK, sizeof(int));

  sprintf(fshape,"%svoid_Shapefinders_iso%f_%d",fname, rho_th1,count_th);
  fp2=fopen(fshape, "w");

  sprintf(fbubble,"%svoid_distribution_iso%f",fname, rho_th1);
  fp3=fopen(fbubble, "w");

  NC = DriveClusterFinder(cluster, rho_th, &ff);

  fclose(fp2); fclose(fp3);

  printf("NC = %d; percolation_no=%d; count_max=%d,  NN_max=%d, count_max_vol=%lf:: vol_max=%lf, NN_max_vol=%d, vol_max_count=%d \n", NC,percolation_no,count_max,NN_max,count_max_vol,vol_max,NN_max_vol,vol_max_count);
  printf("**********************************************\n\n");

  if(NN_check>0){
	  fp = fopen("test_cluster.dat", "w");
	  for(ii=0; ii<NI; ii++){
	    for(jj=0; jj<NJ; jj++) {
	      for(kk=0; kk<NK; kk++)   {
		index = (ii*NJ + jj)*NK + kk;
		//fprintf(fp, "%d\n", cluster[index]);
		if(cluster[index]==NN_check)
		{
			fprintf(fp, "%d\t%d\t%d \t %d \n", ii,jj,kk, cluster[index]);
			//fprintf(fp, "%d\t%d\t%d %d\n", ii,jj,kk,cluster[index]);	
		}
		if(cluster[index]<0){printf("Problem cluster [%d]=%d \n ",index, cluster[index] );} // should not remain -1 anywhere
	      }
	      //fprintf(fp, "\n");
	    }
	   }
	  fclose(fp);
   }


	double ff1=((double)ncpoints/(double) NIJK); double ff_vol=tot_vol/((double) NIJK);
	double LCS=((double) count_max)/((double) count1); double LCS_vol=vol_max/((double) count1); //tot_vol is not used since it depends on the count_th  
  
  printf("count1=%d, ncpoints=%d : they should match. NIJK=%d, ff=%lf, ff1=%lf, ff_vol= %lf \n",count1, ncpoints,NIJK,ff,((double)ncpoints/(double) NIJK),tot_vol/((double) NIJK));

  printf("rho_th1=%lf, iso=%lf \n", rho_th1, iso);
  printf("bad_NE=%d \n",bad_NE);

  //========================================================================================================
  //required in very special cases where we want to study percolation with time evolution.
  //if you don't require this, please comment it out
  sprintf(fcluster_stat,"%sCluster_stat_iso%f",fname, rho_th1);
  fp=fopen(fcluster_stat, "w");
  fprintf(fp,"#NI = %d \t NJ = %d \t NK = %d \n",NI,NJ,NK);
  fprintf(fp,"#rho_th \t NC \t total_count \t NIJK \t FF \t percolation_no \t count_max \t LCS \t NN_max \t count_max_vol \t vol_max \t NN_max_vol \t vol_max_count \n");
  fprintf(fp,"%.16lf \t %d \t %d \t %d \t %.16lf \t %d \t %d \t %.16lf \t %d \t %.16lf \t %.16lf \t %d \t %d \n",rho_th1,NC,count1,NIJK,ff1,percolation_no,count_max,LCS,NN_max,count_max_vol,vol_max,NN_max_vol,vol_max_count);
  fprintf(fp,"\n #bad_NE= %d\n",bad_NE);
  fclose(fp);
  //========================================================================================================

  
  return NC; 
}

//============================ main ends ===================================================================================================
/*********************************** push a element v in the stack S **********************************************************************/
void push(struct Stack **start, int v[5], int x) //x=1: i direction. x=2: j direction . x=3: k direction.
{
	
	if(v[4]>-1){par[x-1]=-1;ijk_shift[x-1]=0;/*printf("...entered special PUSH...dir=%d\n",v[4]); if(v[4]!=x-1){printf("...................................................... PROBLEM....\n \n \n \n \n \n \n \n \n \n \n \n");}*/}
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
		if((ijk[abs(n)-1]+abs(n)/n)>=NIJK_array[abs(n)-1] && nb[abs(n)]==1 && par[abs(n)-1]==1){nb_dummy[4]=abs(n)-1; /*printf("...entered special case...dir=%d\n",abs(n)-1);*/}
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

	if(rho[index]>rho_th) {
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
	/*if(NN==10967){printf("\n.....NN=%d begins......\n\n",NN);}*/
  int nb[4]={ii,1,1,1};par[0]=1;par[1]=1;par[2]=1;//nb[0] stores index. nb[1,2,3]=xp,yp,zp store parity
  int nb_dummy[5]={ii,1,1,1,-1}; int nxt; //nb_dummy[4], i.e. 5th element stores the info to change ijk_shift[required element]=0 in the special case

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

  fprintf(fp3,"%d \t %d \n",NN,count);	
 	
 if(count>=count_th)
 {
 	
 	printf("NC=%d: count= %d :",NN,count); 
  /*if(check_percolation(xmin[NN], xmax[NN], ymin[NN], ymax[NN], zmin[NN], zmax[NN])>0)
  {
  	printf("Percolation happend. NC= %d \n", NN);
  }*/
  
  printf("NC= %d:: ijk_shift= {%d, %d, %d}.....",NN,ijk_shift[0],ijk_shift[1],ijk_shift[2]);
  
  ijk_shift[0]=(ijk_shift[0]+1)%NI; ijk_shift[1]=(ijk_shift[1]+1)%NJ; ijk_shift[2]=(ijk_shift[2]+1)%NK;
  
  printf("NC= %d:: ijk_shift= {%d, %d, %d}...stroing points... \n",NN,ijk_shift[0],ijk_shift[1],ijk_shift[2]);
  
  
   char hname[100]; //char hname1[100];
   FILE *fp; //FILE * fp1;
   if(NN==NN_check){ sprintf(hname, "cluster_%d.dat", NN); //sprintf(hname1, "data_%d.dat", NN);
   fp=fopen(hname, "w");  //fp1=fopen(hname1, "w");
   }
   
  //int * arranged_cluster;
  //arranged_cluster = (int *) calloc(NIJK1, sizeof(int));
  
  int in,jn,kn,i,j,k, index_old, index_new,index_new1; bool clause =true;
  //fprintf(fp1, "#%d\t%d\t%d \t %lf \t %lf \n", NI+1,NJ+1,NK+1, -1215.0,2000.0);
  
  //need the final box of size (NI+1)*(NJ+1)*(NK+1) due to the convention i \equiv i+NI (rajesh sends data in this convention), not \equiv i+NI-1.
  int indi;
  Percolation_happend: indi=0;
  for(in=0;in<_size_x;in++)
  {
  	for(jn=0;jn<_size_y;jn++)
  	{
  		for(kn=0;kn<_size_z;kn++)
  		{
  			i=(in+ijk_shift[0])%NI; j=(jn+ijk_shift[1])%NJ; k=(kn+ijk_shift[2])%NK; //shifting according to ijk_shift[]
  			index_old=get_index(i,j,k); index_new=get_index(in,jn,kn);
  			//arranged_cluster[index_new]=cluster[index_old];
  			index_new1=(in*_size_y+jn)*_size_z+kn;
  			
  			if(cluster[index_old]==NN)
  			{
  				
  				if(NN==NN_check){fprintf(fp, "%d\t%d\t%d \t %d \n", in,jn,kn, cluster[index_old]);}
  				//fprintf(fp1, "%lf \n", rho[index_old]);
  				data[indi]=rho[index_old];
  				if((in==0 || jn==0 || kn==0))
  				{
  					//if the final box faces go through any data point: percolation
  					if(clause==true){
  							printf("Percolation happend.................................NC=%d \n",NN); 
  							clause=false;
  							percolation_no++;
  							//when percolation detected..No shift. Fill data[] from the beginning
  							ijk_shift[0]=0; ijk_shift[1]=0;ijk_shift[2]=0;
  							goto Percolation_happend;
  					} 
  					data[indi]=0.9*rho_th; rho[index_old]=0.9*rho_th;
  				}
  			}
  			else if (cluster[index_old]!=0){data[indi]=0.9*rho_th;}//{fprintf(fp1, "%lf \n", 0.1*rho_th);}
  			else {data[indi]=rho[index_old];}//{fprintf(fp1, "%lf \n", rho[index_old]);}
  			indi++;
  			
  		}
  	}
  }
  
 
   
   if(NN==NN_check){fclose(fp);}//fclose(fp1);
   
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
    		
    		if(L>=B && B>=T){printf("L>=B>=T: ");} 
    		
    		else 
    		{
    			printf("L,B,T not in order: ");
    			//if(Shape3<Shape2){L=Shape2; B=Shape3;}
    			
    		}
    		double P=(B-T)/(B+T); double F=(L-B)/(L+B); double F1=(L1-B)/(L1+B);
    
    		printf(" NN=%d, count=%d, Vol=%lf, area=%lf, imc=%lf, imc1=%lf, imc2=%lf, genus=%lf, T=%lf, B=%lf, L=%lf, Planarity=%lf, Filamentarity=%lf, Filamentarity1=%lf \n",NN,count,vol,area,imc,imc1,imc2,genus,T,B,L,P,F,F1);
  		printf("*********************************************************************************************************************\n\n");
  		
  		fprintf(fp2,"%d \t %d \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \n", NN,count,vol,area,imc,imc1,imc2,genus,T,B,L,P,F,F1);
  		tot_vol=tot_vol+vol; //tot_vol only includes count>some value set
  		if(vol>vol_max){vol_max=vol; NN_max_vol=NN; vol_max_count=count;}
  		if (count==count_max){count_max_vol=vol;}
  		//clean_all();
  		
  		
  	}
  return NN+1;
}
