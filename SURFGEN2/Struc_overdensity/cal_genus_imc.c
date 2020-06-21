// Needs fixup.c. Used ConnectedVertexNode tree to store each connected vertex to each VertexNode.

# include "MarchingCube.h"
#include "imc.h"
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

/*//char fin[100]="test_tri.d";
char fin[100]="../SHAPEFINDER_pbc/tri_1.dat";
//char fin[100]="../New_TEST/tri_1.dat";*/
int extVert;
int NV; 
int max_vert_shared=100;



double magnitude2(Vector V)
{
	double mag2=V.x*V.x+V.y*V.y+V.z*V.z;
	return mag2;
	//return sqrt(mag2);
}


Vector add_Vec(Vector V1, Vector V2)
{
	Vector V;
	V.x=V1.x+V2.x; V.y=V1.y+V2.y; V.z=V1.z+V2.z;
	return V;
}

Vector sub_Vec(Vector V1, Vector V2) //gives you V1-V2
{
	Vector V;
	V.x=V1.x-V2.x; V.y=V1.y-V2.y; V.z=V1.z-V2.z;
	return V;
}

Vector null_Vec()
{
	Vector V;
	V.x=0.0; V.y=0.0; V.z=0.0;
	return V;
}

//=============================== structure ends ======================================================================================
/*
int no_of_lines() //for C
{
	FILE *fp;
    	int count = 0; char c;
    	fp = fopen(fin, "r");
    	if (fp == NULL)
    	{
        	printf("Could not open file %s", fin);
        	return 0;
    	}
 
    	// Extract characters from file and store in character c
    	for (c = getc(fp); c != EOF; c = getc(fp))
        	if (c == '\n') // Increment count if this character is newline
            		count = count + 1;
 
    		// Close the file
    		fclose(fp);
    		return count;



}

void declare_array()
{
        
        TV = (double ***)malloc(NT*sizeof(double**));
        TNV= (Vector **)malloc(NT*sizeof(Vector*));
        int i;
        for ( i = 0; i<= NT; i++) // to avoid segmentation fault, <=NT is used. Actually <NT creates an array of size NT. Still I am getting the segmentation fault in initialize(), in last of while. This problem does not appear if we use an array TV[NT][3][3] instead!!
        {
                TV[i] = (double **) malloc(3*sizeof(double *));
                TNV[i] = (Vector *) malloc(3*sizeof(Vector)); 
                int j;
                for ( j = 0; j < 3; j++) 
                {
                        TV[i][j] = (double *)malloc(3*sizeof(double));
                }

       }
}

void initialize()
{
        int i; int ri=0;
        FILE *fp;
    	fp = fopen(fin, "r");
        
        i = fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &TV[ri][0][0], &TV[ri][0][1], &TV[ri][0][2], &TV[ri][1][0], &TV[ri][1][1], &TV[ri][1][2], &TV[ri][2][0], &TV[ri][2][1], &TV[ri][2][2]); ri++;
        while (i != EOF)
        {
           	i=fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &TV[ri][0][0], &TV[ri][0][1], &TV[ri][0][2], &TV[ri][1][0], &TV[ri][1][1], &TV[ri][1][2], &TV[ri][2][0], &TV[ri][2][1], &TV[ri][2][2]);
           	if(i != EOF){ri++;}
        }

        fclose(fp);

	//printf("ri= %d \n",ri);

    
}*/
//=========================================================================================================================================


//******************************* tree algorithm starts *********************************************************************************** 

Vector get_normal(int tn)
{
	Vector Ver1, Ver2,nor;
	double vx=TV[tn][0][0]; double vy=TV[tn][0][1]; double vz=TV[tn][0][2];
	Ver1.x=TV[tn][1][0]-vx; Ver1.y=TV[tn][1][1]-vy; Ver1.z=TV[tn][1][2]-vz;
	Ver2.x=TV[tn][2][0]-vx; Ver2.y=TV[tn][2][1]-vy; Ver2.z=TV[tn][2][2]-vz;
	nor=CROSSPRODUCT(Ver2,Ver1); //changed order to arrange the normals outward of the surface, compitable with Prakashda's storing
	nor=Normalise(nor);
	return nor;
}
 
double get_ep_theta(int tn1, int tn2)//(int ver_no1, int ver_no2)
{
	//int tn1=ver_no1/3; int vn1=ver_no1%3;
	//int tn2=ver_no2/3; int vn2=ver_no2%3;
	Vector nor1, nor2,c2c;
	//Vector t1,t2;
	double Px1,Py1,Pz1,Px2,Py2,Pz2; //stores coordinates (x,y,z) of centres of the two triangles
	nor1=get_normal(tn1); nor2=get_normal(tn2);
	double costheta=DOTPRODUCT(nor1,nor2);
	if(fabs(costheta)>1.000001){printf("BAD...dot_product=%lf\n",costheta);}
	if(costheta>1.0){costheta=1.0;}
	if(costheta<-1.0){costheta=-1.0;}
	double theta=acos(costheta); double ep;
	//printf("theta=%lf \n",theta);
	//if(theta!=theta){printf("BAD....theta=%lf dot_product=%lf\n",theta,costheta);}
	double I1,I2;
	
//**************** for first triangle tn1 ********************************************	
	
	Px1=(TV[tn1][0][0]+TV[tn1][1][0]+TV[tn1][2][0])/3.0; Py1=(TV[tn1][0][1]+TV[tn1][1][1]+TV[tn1][2][1])/3.0; Pz1=(TV[tn1][0][2]+TV[tn1][1][2]+TV[tn1][2][2])/3.0;
	
	Px2=(TV[tn2][0][0]+TV[tn2][1][0]+TV[tn2][2][0])/3.0; Py2=(TV[tn2][0][1]+TV[tn2][1][1]+TV[tn2][2][1])/3.0; Pz2=(TV[tn2][0][2]+TV[tn2][1][2]+TV[tn2][2][2])/3.0;
	
	c2c.x=Px2-Px1; c2c.y=Py2-Py1; c2c.z=Pz2-Pz1; //vector fron center of tn1 to centre of tn2
	
	/*t1.x=Px1-x1; t1.y=Py1-y1; t1.z=Pz1-z1;
	if(fabs(DOTPRODUCT(nor1,t1))>0.000001){printf("PROBLEM.... P.nor1=%lf \n",DOTPRODUCT(nor1,t1));}*/
	
	
	I1=DOTPRODUCT(nor1,c2c);
	
//**************** for second triangle tn2 ********************************************
	
	c2c.x=-c2c.x;c2c.y=-c2c.y;c2c.z=-c2c.z; //vector fron center of tn2 to centre of tn1
	
	/*t2.x=Px2-x1; t2.y=Py2-y1; t2.z=Pz2-z1;
	if(fabs(DOTPRODUCT(nor2,t2))>0.000001){printf("PROBLEM.... P.nor2=%lf \n",DOTPRODUCT(nor2,t2));}*/
	
	//P1.x=nor1.x+TV[tn1][0][0]; P1.y=nor1.y+TV[tn1][0][1]; P1.z=nor1.z+TV[tn1][0][2];
	I2=DOTPRODUCT(nor2,c2c);
//*************************************************************************************	
	if(I1>0 && I2>0){ep=-1.0; /*printf("Concave..\n");*/}
	else if (I1<0 && I2<0){ep=1.0;/*printf("..Convex\n");*/}
	else
	{
		//should not enter here: Can enter only when I1,I2 are very close to zero
		//printf("BAD............. I1=%lf, I2=%lf--tn=%d vn=%d, tree node: tnc=%d, vnc=%d\n",I1, I2,tn1,vn1,tn2,vn2);
		if(fabs(I1*I2)>0.00001){printf("BAD............. I1=%lf, I2=%lf--tn=%d, tree node: tnc=%d \n",I1, I2,tn1,tn2);}
		ep=0.0; //can be assigned ep=0: because theta should be very close to zero too;
	}
	return theta*ep;
} 
  
 void make_edge2(int tn, int vd, struct VertexNode * ptr)
 {
	int ver_no=3*tn+vd;
 	//if(ptr->root_cvert!=NULL){printf("Not Null...\n");}
  	insert_conn_vert(ptr,ver_no); //calculation of NE (NE1) happens in insert_conn_vert() and createConnectedVertNode()
  	//printf("poniter...\n");
 
 } 
  
  void make_edge1(int tn, int vn, struct VertexNode * ptr) //calls void make_edge2
  {
  	int va,vb;
  	if(vn==0){va=1;vb=2;}
  	else if(vn==1){va=2;vb=0;}
  	else{va=0; vb=1;}	
  	make_edge2(tn,va,ptr);
  	make_edge2(tn,vb,ptr);
  		
  }

  struct VertexNode * createNode(int tn, int vn, Vector nvw) 
  {
  	NV++;
  	double x,y,z;
  	x=TV[tn][vn][0];y=TV[tn][vn][1];z=TV[tn][vn][2];
        struct VertexNode *newnode;
        newnode = (struct VertexNode *)malloc(sizeof(struct VertexNode));
        newnode->x = x;
        newnode->y = y;
        newnode->z = z;
        newnode->root_cvert=NULL;
        newnode->conn_tri_no=1;
        newnode->nv=nvw;
        int ver_no=3*tn+vn;
        newnode-> shared_tri[0]=ver_no;
        TNV[tn][vn]=nvw;
        
        make_edge1(tn,vn,newnode);
        
        newnode->color = RED;
        newnode->link[0] = newnode->link[1] = NULL;
        return newnode;
  }
  
 //void fixup_RBT(struct VertexNode **, int*, int); //called by insertion
//*********************** insert in the tree ******************************************************************************************
  void insertion (int tn, int vn, Vector nvw) 
  {
  	//int data;
  	double x,y,z;
  	Vector new_nv;
  	x=TV[tn][vn][0];y=TV[tn][vn][1];z=TV[tn][vn][2];
  	
        struct VertexNode *stack[98], *ptr, *newnode;
        int dir[98], ht = 0, index;
        ptr = root;
        if (!root) 
        {
                root = createNode(tn,vn,nvw);
                return;
        }
        stack[ht] = root;
        dir[ht++] = 0;
        /* find the place to insert the new node */
        while (ptr != NULL) {
                if (ptr->x == x && ptr->y == y && ptr->z == z) 
                {
                        extVert++;// printf("Duplicates Not Allowed!!\n");
                        
                        make_edge1(tn,vn,ptr);
                        
                        new_nv=add_Vec(ptr->nv,nvw);
                        ptr->nv=new_nv;
                        int ver_no=3*tn+vn; ptr-> shared_tri[ptr->conn_tri_no]=ver_no;
                        
                        ptr->conn_tri_no=ptr->conn_tri_no+1;
                        int i=0;
                        for(i=0;i<ptr->conn_tri_no;i++)
                        {
                        	int tnc=ptr->shared_tri[i]/3; int vnc=ptr->shared_tri[i]%3;
                        	TNV[tnc][vnc]=new_nv;
                        }
                        
                        return;
                }
                //index = (data - ptr->data) > 0 ? 1 : 0;
                if(ptr->x>x) {index=1;}
                else if (ptr->x<x){index=0;}
                else
                {
                	if(ptr->y>y){index=1;}
                	else if(ptr->y<y){index=0;}
                	else
                	{
                		if(ptr->z>z){index=1;}
                		else if(ptr->z<z){index=0;}
                		else
                		{printf("\n This is BAD. Should not come \n \n");
                		
                			//extVert++;// printf("Duplicates Not Allowed!!\n");
                        		//return;
                		}
                	
                	}
                }
                
                
                stack[ht] = ptr;
                ptr = ptr->link[index];
                dir[ht++] = index;
        }
        /* insert the new node */
        stack[ht - 1]->link[index] = newnode = createNode( tn,vn,nvw); 
        fixup_VertexNode(stack,dir,ht);
        
}

//*********************************************************************************************************************************
//#################################################################################################################################

Vector tri_nor_weight(int tn, int vn)
{
	Vector nv, nor,Ver1,Ver2; int v1,v2;
	if(vn==0){v1=1;v2=2;} 
	if(vn==1){v1=2;v2=0;} 
	if(vn==2){v1=0;v2=1;}
	double vx=TV[tn][vn][0]; double vy=TV[tn][vn][1]; double vz=TV[tn][vn][2];
	Ver1.x=TV[tn][v1][0]-vx; Ver1.y=TV[tn][v1][1]-vy; Ver1.z=TV[tn][v1][2]-vz;
	Ver2.x=TV[tn][v2][0]-vx; Ver2.y=TV[tn][v2][1]-vy; Ver2.z=TV[tn][v2][2]-vz;

	nor=CROSSPRODUCT(Ver2,Ver1); //changed order to arrange the normals outward of the surface, compitable with Prakashda's storing
	double mag2Ver1=magnitude2(Ver1); double mag2Ver2=magnitude2(Ver2);
	double c=mag2Ver1*mag2Ver2;
	
	nv.x=nor.x/c; nv.y=nor.y/c; nv.z=nor.z/c;
	return nv;
	
}

void print_TNV()
{
	int i,j;
	for(i=0;i<NT;i++)
	{
		printf("[");
		for(j=0;j<3;j++)
		{
			printf("( %lf, %lf, %lf),",TNV[i][j].x,TNV[i][j].y,TNV[i][j].z);
		}
		printf("]\n");
	}

}

void print_Vec(Vector V)
{
	printf ("(%lf, %lf, %lf)  ",V.x,V.y,V.z);
}

double inv_MTM[3][3];

void inverse_MTM(double MTM[3][3])
{
	double detMTM=MTM[0][0]*(MTM[1][1]*MTM[2][2]-MTM[1][2]*MTM[2][1])+MTM[0][1]*(MTM[1][2]*MTM[2][0]-MTM[1][0]*MTM[2][2])+MTM[0][2]*(MTM[1][0]*MTM[2][1]-MTM[1][1]*MTM[2][0]);
	
	if(detMTM==0 || detMTM!=detMTM){printf("Problem:: detMTM=%lf \n",detMTM);}
	
	inv_MTM[0][0]=MTM[1][1]*MTM[2][2]-MTM[1][2]*MTM[2][1]; inv_MTM[0][0]=inv_MTM[0][0]/detMTM;
	inv_MTM[0][1]=MTM[0][2]*MTM[2][1]-MTM[0][1]*MTM[2][2]; inv_MTM[0][1]=inv_MTM[0][1]/detMTM;
	inv_MTM[0][2]=MTM[0][1]*MTM[1][2]-MTM[0][2]*MTM[1][1]; inv_MTM[0][2]=inv_MTM[0][2]/detMTM;
	inv_MTM[1][0]=MTM[1][2]*MTM[2][0]-MTM[1][0]*MTM[2][2]; inv_MTM[1][0]=inv_MTM[1][0]/detMTM;
	inv_MTM[1][1]=MTM[0][0]*MTM[2][2]-MTM[0][2]*MTM[2][0]; inv_MTM[1][1]=inv_MTM[1][1]/detMTM;
	inv_MTM[1][2]=MTM[0][2]*MTM[1][0]-MTM[0][0]*MTM[1][2]; inv_MTM[1][2]=inv_MTM[1][2]/detMTM;
	inv_MTM[2][0]=MTM[1][0]*MTM[2][1]-MTM[1][1]*MTM[2][0]; inv_MTM[2][0]=inv_MTM[2][0]/detMTM;
	inv_MTM[2][1]=MTM[0][1]*MTM[2][0]-MTM[0][0]*MTM[2][1]; inv_MTM[2][1]=inv_MTM[2][1]/detMTM;
	inv_MTM[2][2]=MTM[0][0]*MTM[1][1]-MTM[0][1]*MTM[1][0]; inv_MTM[2][2]=inv_MTM[2][2]/detMTM;
	
	/*int i,j,k; double sum;
	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			sum=0;
			for(k=0;k<6;k++)
			{
				sum=sum+MTM[i][k]*inv_MTM[k][j];	
			}
			printf("%lf ",sum);
		}	
		printf("\n");
	}*/
	
}



double mean_curv(int tn)
{
	//printf("\n\n");
	Vector e0, e1, e2;
	Vector dn10, dn21, dn02;
	
	dn10=sub_Vec(TNV[tn][1],TNV[tn][0]); dn21=sub_Vec(TNV[tn][2],TNV[tn][1]); dn02=sub_Vec(TNV[tn][0],TNV[tn][2]);
	e0.x=TV[tn][2][0]-TV[tn][1][0]; e0.y=TV[tn][2][1]-TV[tn][1][1]; e0.z=TV[tn][2][2]-TV[tn][1][2];
	e1.x=TV[tn][0][0]-TV[tn][2][0]; e1.y=TV[tn][0][1]-TV[tn][2][1]; e1.z=TV[tn][0][2]-TV[tn][2][2];
	e2.x=TV[tn][1][0]-TV[tn][0][0]; e2.y=TV[tn][1][1]-TV[tn][0][1]; e2.z=TV[tn][1][2]-TV[tn][0][2];
	
	//define (u,v), the coordinate system on the face plane in terms of my original cordinates (i,j,k).
	Vector u,v,v_temp,N;
	u=Normalise(e0);
	N=CROSSPRODUCT(e0,e1);
	v_temp=CROSSPRODUCT(N,e0);
	v=Normalise(v_temp);
	//printf("%lf \n",DOTPRODUCT(u,v));
	
	
	double x1,x2,x3,x4,x5,x6,y1,y2,y3,y4,y5,y6;
	double imc_tn;
	x1=DOTPRODUCT(e0,u); x2=DOTPRODUCT(e0,v); x3=DOTPRODUCT(e1,u);x4=DOTPRODUCT(e1,v); x5=DOTPRODUCT(e2,u); x6=DOTPRODUCT(e2,v);
	y1=DOTPRODUCT(dn21,u); y2=DOTPRODUCT(dn21,v); y3=DOTPRODUCT(dn02,u);y4=DOTPRODUCT(dn02,v); y5=DOTPRODUCT(dn10,u); y6=DOTPRODUCT(dn10,v);
	
	double M[6][3]; double MT[3][6]; double MTM[3][3]; 
	M[0][0]=MT[0][0]=x1; M[0][1]=MT[1][0]=0.; M[0][2]=MT[2][0]=x2;
	M[1][0]=MT[0][1]=0.; M[1][1]=MT[1][1]=x2; M[1][2]=MT[2][1]=x1;
	M[2][0]=MT[0][2]=x3; M[2][1]=MT[1][2]=0.; M[2][2]=MT[2][2]=x4;
	M[3][0]=MT[0][3]=0.; M[3][1]=MT[1][3]=x4; M[3][2]=MT[2][3]=x3;
	M[4][0]=MT[0][4]=x5; M[4][1]=MT[1][4]=0.; M[4][2]=MT[2][4]=x6;
	M[5][0]=MT[0][5]=0.; M[5][1]=MT[1][5]=x6; M[5][2]=MT[2][5]=x5;
	
	int i,j,k; double sum;
	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			sum=0;
			for(k=0;k<6;k++)
			{
				sum=sum+MT[i][k]*M[k][j];
			}
			MTM[i][j]=sum;
		}	
	}
	
	inverse_MTM(MTM); //inverse is stored in gobal array inv_MTM[3][3]
	
	double X[3], Y[6]; Y[0]=y1; Y[1]=y2; Y[2]=y3; Y[3]=y4;Y[4]=y5; Y[5]=y6; X[0]=X[1]=X[2]=0;
	
	
	for(i=0;i<3;i++)
	{
		
		for(j=0;j<3;j++)
		{
			sum=0;
			for(k=0;k<6;k++)
			{
				sum=sum+MT[j][k]*Y[k];
			}
			X[i]=X[i]+inv_MTM[i][j]*sum;
		}
			
	}
	double trace=X[0]+X[1];
	double area=0.5*sqrt(magnitude2(N));
	imc_tn=0.5*area*trace; 
	return imc_tn;
	
	/*
	//printf(" tn= %d %lf %lf %lf\n",tn,magnitude2(e0),magnitude2(dn21),magnitude2(v));
	if(tn==1 || tn==3){printf("tn =%d  ",tn); print_Vec(e0); print_Vec(dn21);print_Vec(v); printf("\n \n ");}
	double temp=x2*x3-x1*x4; 
	A=(y3*x2-y1*x4)/temp;
	C=(y1*x3-y3*x1)/temp;
	B=(y2*x2*x3-y2*x1*x4-y1*x3*x1+y3*x1*x1)/(x2*temp);
	
	//if(x2==0){printf("x2=0. tn= %d \n",tn);}
	//printf("%lf %lf %lf %lf\n",x1,x2,x3,x4);
	//printf("%lf %lf %lf \n",y1,y2,y3);
	//printf("%lf %lf %lf \n",A,B,C);
	
	double T=A*B; double D=A*B-C*C;
	double check=T*T-4.0*D;
	//if(check<0){printf("Problem: T*T-4.0*D<0: tn= %d",tn);}
	
	double k1=0.5*T +0.5*sqrt(check);
	double k2=0.5*T -0.5*sqrt(check);
	double area=0.5*sqrt(magnitude2(N));
	imc_tn=0.5*area*T;
	//printf("tn= %d, imc_tn= %lf, k1=%lf, k2=%lf, check=%lf \n",tn, imc_tn,k1,k2,check);
	return imc_tn;*/
	
}

double max_length(int tn)
{
	double x1,y1,z1,x2,y2,z2,x3,y3,z3;
	
	x1=TV[tn][0][0];y1=TV[tn][0][1];z1=TV[tn][0][2];x2=TV[tn][1][0];y2=TV[tn][1][1];z2=TV[tn][1][2];x3=TV[tn][2][0];y3=TV[tn][2][1];z3=TV[tn][2][2];
	
	double len, len_maxi=0;
	len=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1));
	if(len>len_maxi){len_maxi=len;}
	len=sqrt((x2-x3)*(x2-x3)+(y2-y3)*(y2-y3)+(z2-z3)*(z2-z3));
	if(len>len_maxi){len_maxi=len;}
	len=sqrt((x1-x3)*(x1-x3)+(y1-y3)*(y1-y3)+(z1-z3)*(z1-z3));
	if(len>len_maxi){len_maxi=len;}
	return len_maxi;
	
}

double compute_genus_imc()
{
  
	root=NULL; NE=0; NE1=0; imc1=0;imc2=0;
	NV=0; extVert=0;  	
  
  	/*NT=no_of_lines();
        declare_array(); 
	initialize();
	printf("NT=%d",NT);*/
        
        Vector nvw;
        int tn=0; int vn=0;
        //printf("************* for loop starts ************************************** \n");
        for(tn=0;tn<NT;tn++)
        {
        	for(vn=0;vn<3;vn++)
        	{
        		nvw=tri_nor_weight(tn,vn);
        		//x=TV[tn][vn][0];y=TV[tn][vn][1];z=TV[tn][vn][2];
        		insertion(tn,vn,nvw);
        	}
        	
        }
        
        
        NE=NE/2.0; NE1=NE1/2.0;
        //double max_len;
        double genus= 1.0-0.5*(NT-NE+NV);
        if(NE!=1.5*NT) {bad_NE++;}
        printf(" \n NT= %d, NV= %d, extVert= %d, NE=%lf, NE1=%lf, diff_NE=%lf, genus= %lf \n",NT,NV,extVert,NE,NE1,1.5*NT-NE,genus);

        for(tn=0;tn<NT;tn++)
        {
        	for(vn=0;vn<3;vn++)
        	{
        		TNV[tn][vn]=Normalise(TNV[tn][vn]);
        	}
        	
        	imc1=imc1+mean_curv(tn);
        	//max_len=max_length(tn);
        	//if(max_len>len_max2){len_max2=max_len;}
        }
        //print_TNV();
        imc2=imc2/2.0;
        //printf("imc1= %lf, imc2=%lf \n", imc1,imc2);
        
        free_VertexNode(root);
        root=NULL;
        //printf("reached root cleaning in cal_genus_imc.c \n");
        return genus;
        //return 0;
  }
