# include "MarchingCube.h"
# include "LookUpTable.h"
# define N 10000
# define MAXLENGTH 15
# define ALLOC_SIZE 65536  //what is so special about 65536 ??
# define EPSILON 1.e-05
# define BYTE2MB 1./(1024.*1024.)
# define NEXT(i) ((i)<2 ? (i)+1 : (i)-2)
# define PREV(i) ((i)>0 ? (i)-1 : (i)+2)

/*defining Global variable*/

bool _originalMC = FALSE;
int _Ntrigs = 0, _Nverts = 0;
int _nverts = 0, _ntrigs = 0;
int _case, _config, _subconfig;
int _i, _j, _k, _lut_entry;
int *_x_verts = NULL, *_y_verts = NULL, *_z_verts = NULL;
/* int *_xmin = NULL, *_xmax = NULL, *_ymin = NULL, *_ymax = NULL, *_zmin = NULL, *_zmax = NULL; */
//int *_clusters=NULL;
//double *data = NULL, *rho = NULL;
double _cube[8];
Triangle *_triangles=NULL;
Vertex *_vertices=NULL;
int tunnelOrientation = 0;
Vector *Normals=NULL;
Vector *pdir1=NULL, *pdir2=NULL;
int unused __attribute__((unused));
double Niso;

int count = 0;

//int NT2=0; //for genus calculation; NE and NV are in cal_genus.c //NT1 is calculated in store_triangle() in Merchincube.c
//we name it NT1 because NT is also used globaly in some functions


/*---------------------debug tools-----------------------------*/
void print_cube() {
  printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", _cube[0], _cube[1],
         _cube[2], _cube[3], _cube[4], _cube[5], _cube[6], _cube[7]);
}
void print_vertex(int ii) {
  printf("%e\t%e\t%e\n", _vertices[ii].x, _vertices[ii].y, _vertices[ii].z);
}

void printsize() {
  printf("************ total ram used upto this******** \n ");
  printf("vertices(%f MB), Triangle(%f MB) \n", malloc_usable_size(_vertices)/(1024.*1024.), malloc_usable_size(_triangles)/(1024.*1024.));
  printf("clusters(%f MB), data(%f MB), rho(%f MB)\n", malloc_usable_size(_clusters)/(1024.*1024.), malloc_usable_size(data)/(1024.*1024.), malloc_usable_size(data)/(1024.*1024.));
  printf("************ total ram used upto this******** \n \n");

}

/*--------------------Initialization---------------------------*/
void init_temps() { //used in init_all
  int ii, jj, kk, index;

  _x_verts=(int*)calloc(_size_x*_size_y*_size_z, sizeof(int));
  _y_verts=(int*)calloc(_size_x*_size_y*_size_z, sizeof(int));
  _z_verts=(int*)calloc(_size_x*_size_y*_size_z, sizeof(int));

  for(ii=0; ii<_size_x; ii++)
    for(jj=0; jj<_size_y; jj++)
      for(kk=0; kk<_size_z; kk++) {
        index = kk + (jj * _size_z) + (ii * _size_z * _size_y);
        _x_verts[index]=-1;_y_verts[index]=-1;_z_verts[index]=-1;
      }
}

void init_all()  {   //used in PrepareGrid function
  init_temps() ;

  _nverts = _ntrigs = 0 ;
  _Nverts = _Ntrigs = ALLOC_SIZE ;  //what is signifiacs of ALLOC_SIZE?

  _vertices=(Vertex *)calloc(_Nverts, sizeof(Vertex));  //Vertex structure variable has coordinates (x,y,z) and normals (nx,ny,nz)
  _triangles=(Triangle *)calloc(_Ntrigs, sizeof(Triangle)); //Triangle structure variable has int v1,v2,v3
}

void clean_temps()  { //only called in clean_all
  free(_x_verts);
  free(_y_verts);
  free(_z_verts);

  _x_verts = NULL;
  _y_verts = NULL;
  _z_verts = NULL;
}

void clean_all()  {  //used in main.c at the end.
  clean_temps();
  free(_vertices);
  free(_triangles);
  _nverts = _ntrigs = 0 ;
  _Nverts = _Ntrigs = 0 ;

 // _size_x = _size_y = _size_z = -1 ;
}

void free_TV_TNV()
{
	int i, j;
	for(i=0;i<NT;i++)
	{
		for(j=0;j<3;j++)
		{
			free(TV[i][j]);
		}
		free(TV[i]);
		free(TNV[i]);
	}

	free(TV);
	//TV=NULL;

	free(TNV);
	//TNV=NULL;	

}

/*------------------------main algorithm-------------------------*/
void run(int NN) {  //run is used in main.c once
  int p;
  
  init_all();
	
  compute_intersection_points();

  for( _i = 0 ; _i < _size_x-1 ; _i++ )
    for( _j = 0 ; _j < _size_y-1 ; _j++ )
      for( _k = 0 ; _k < _size_z-1 ; _k++ )  {
        _lut_entry = 0;
        for(p = 0; p < 8; ++p)    {
          _cube[p] = get_data(_i+((p^(p>>1))&1), _j+((p>>1)&1),
                              _k+((p>>2)&1)) - iso;
          if(fabs( _cube[p]) < DBL_EPSILON) _cube[p] = 0;   //FLT_EPSILON =1.19209 e-7; why is so large?
          if(_cube[p] >=0) _lut_entry += 1 << p;
        }
        if(_lut_entry>0) {
          process_cube();
        }
      }
      
      no_of_tri=store_triangle(NN); 

    	area=compute_area();
 	//printf("area_done..\n");
 	vol=compute_vol();
 	//printf("volume_done..\n");
    	//imc=compute_imc(); //not needed

 		imc=0.0;
    	//imc = compute_curvature();
    	//printf("imc_done..\n");
    	clean_all();
	//printf("cleaning done..\n");
    	genus=compute_genus_imc();   //genus=compute_genus() for cal_genus.c
    	//clean_all();
    	//printf("genus_done..\n");
    	free_TV_TNV();
    	//printf("free TV_TNV done..\n");
}

void compute_intersection_points() {  //called in run only. run() is used in main.c once
  int index, tempi;
  int *flag;

  flag=(int*)calloc(_size_x*_size_y*_size_z, sizeof(int));

  for(tempi=0; tempi<_size_x*_size_y*_size_z; tempi++)
    flag[tempi]=-1; // all elements of flag array has -1


  for(_k = 0; _k < _size_z; _k++)
    for(_j = 0; _j < _size_y; _j++)
      for(_i = 0; _i < _size_x; _i++) {  //what is _cube[8]
        _cube[0] = get_data(_i, _j, _k) - iso; //_cube is double _cube[8]

        if(_i < _size_x - 1)
          _cube[1] = get_data(_i+1, _j, _k) - iso;
        else
          _cube[1] = _cube[0];
        if(_j < _size_y - 1)
          _cube[3] = get_data(_i,_j+1, _k) - iso; //where is cube 2?
        else
          _cube[3] = _cube[0];
        if(_k < _size_z - 1)
          _cube[4] = get_data(_i, _j,_k+1) - iso;
        else
          _cube[4] = _cube[0];
        if(fabs(_cube[0])<DBL_EPSILON) _cube[0]=0.;
        if(fabs(_cube[1])<DBL_EPSILON) _cube[1]=0.;
        if(fabs(_cube[3])<DBL_EPSILON) _cube[3]=0.;
        if(fabs(_cube[4])<DBL_EPSILON) _cube[4]=0.;

        /* define C now to check for intersection point */

        if(_cube[0]>0)     {
          //done for x-direction
          index=_k + _j*_size_z + (_i+1)*_size_z*_size_y;
          if(_cube[1]<0)   {
            set_x_vert(add_x_vertex(), _i, _j, _k) ;
            //printf("1.1 %d\t", flag[index]);
          }
          else if(fabs(_cube[1])<DBL_EPSILON)   {
            if(flag[index]==-1) {
              tempi=add_x_vertex();
              set_x_vert(tempi, _i, _j, _k ) ;
              flag[index] = tempi;
            }
            else
              set_x_vert(flag[index], _i, _j, _k);
            //  printf("1.1.1 %d\t", flag[index]);
          }
          //done for y-direction
          index=_k + (_j+1)*_size_z + _i*_size_z*_size_y;
          if(_cube[3]<0)   {
            set_y_vert(add_y_vertex(), _i, _j, _k) ;
            //printf("1.2 %d\t", flag[index]);
          }
          else if(fabs(_cube[3])<DBL_EPSILON)   {
            if(flag[index]==-1)  {
              tempi=add_y_vertex();
              set_y_vert(tempi, _i, _j, _k ) ;
              flag[index] = tempi;
            }
            else
              set_y_vert(flag[index], _i, _j, _k);
            //printf("1.2.1 %d\t", flag[index]);
          }

          //done for z direction
          index=_k+1 + _j*_size_z + _i*_size_z*_size_y;
          if(_cube[4]<0)   {
            set_z_vert(add_z_vertex(), _i, _j, _k) ;
            //printf("1.3 %d\t", flag[index]);
          }
          else if(fabs(_cube[4])<DBL_EPSILON)  {
            if(flag[index]==-1)  {
              tempi=add_z_vertex();
              set_z_vert(tempi, _i, _j, _k ) ;
              flag[index] = tempi;
            }
            else
              set_z_vert(flag[index], _i, _j, _k);
            //printf("1.3.1 %d\t", flag[index]);
          }

        }
        else if(_cube[0]<0)   {

          //done for x-direction
          index=_k + _j*_size_z + (_i+1)*_size_z*_size_y;
          if(_cube[1]>0)  {
            set_x_vert(add_x_vertex(), _i, _j, _k) ;
            //printf("2.1 %d\t", flag[index]);
          }
          else if(fabs(_cube[1])<DBL_EPSILON)  {
            if(flag[index]==-1)
              {
                tempi=add_x_vertex();
                set_x_vert(tempi, _i, _j, _k ) ;
                flag[index] = tempi;
              }
            else
              set_x_vert(flag[index], _i, _j, _k);
            //printf("2.1.1 %d\t", flag[index]);
          }

          //done for y-direction
          index=_k + (_j+1)*_size_z + _i*_size_z*_size_y;
          if(_cube[3]>0) {
            set_y_vert(add_y_vertex(), _i, _j, _k) ;
            //printf("2.2 %d\t", flag[index]);
          }
          else if(fabs(_cube[3])<DBL_EPSILON)  {
            if(flag[index]==-1)
              {
                tempi=add_y_vertex();
                set_y_vert(tempi, _i, _j, _k ) ;
                flag[index] = tempi;
              }
            else
              set_y_vert(flag[index], _i, _j, _k);
            //printf("2.2.1 %d\t", flag[index]);
          }
          //done for z direction
          index=_k+1 + _j*_size_z + _i*_size_z*_size_y;
          if(_cube[4]>0)  {
            set_z_vert(add_z_vertex(), _i, _j, _k) ;
            //printf("%d\t%d\t%d\t%e\t%e\n", _i, _j, _k, get_data(_i, _j, _k), get_data(_i, _j, _k+1));
            //printf("2.3 %d\t", flag[index]);
          }
          else if(fabs(_cube[4])<DBL_EPSILON)  {
            if(flag[index]==-1)  {
              tempi=add_z_vertex();
              set_z_vert(tempi, _i, _j, _k ) ;
              flag[index] = tempi;
            }
            else
              set_z_vert(flag[index], _i, _j, _k);
            //printf("2.3.1 %d\t", flag[index]);
          }

        }
        else  {
          index=_k + _j*_size_z + _i*_size_z*_size_y;
          if(_cube[1]>0) {
            if(flag[index]!=-1)
              set_x_vert(flag[index], _i, _j, _k);
            else  {
              flag[index]  = add_x_vertex();
              set_x_vert(flag[index], _i, _j, _k);
            }
            //printf("3.1.1 %d\t", flag[index]);
          }
          else if(_cube[1]<0) {
            if(flag[index]!=-1)
              set_x_vert(flag[index], _i, _j, _k);
            else  {
              flag[index]  = add_x_vertex();
              set_x_vert(flag[index], _i, _j, _k);
            }
            //printf("3.1.2 %d\n", flag[index]);
          }

          if(_cube[3]>0) {
            if(flag[index]!=-1)
              set_y_vert(flag[index], _i, _j, _k);
            else  {
              flag[index]  = add_y_vertex();
              set_y_vert(flag[index], _i, _j, _k);
            }
            //printf("3.2.1 %d\t", flag[index]);
          }
          else if(_cube[3]<0) {
            if(flag[index]!=-1)
              set_y_vert(flag[index], _i, _j, _k);
            else  {
              flag[index]  = add_y_vertex();
              set_y_vert(flag[index], _i, _j, _k);
            }
            //printf("3.2.2 %d\t", flag[index]);
          }
          if(_cube[4]>0) {
            if(flag[index]!=-1)
              set_z_vert(flag[index], _i, _j, _k);
            else  {
              flag[index]  = add_z_vertex();
              set_z_vert(flag[index], _i, _j, _k);
            }
            //printf("3.3.1 %d\t", flag[index]);
          }
          else if(_cube[4]<0) {
            if(flag[index]!=-1)
              set_z_vert(flag[index], _i, _j, _k);
            else  {
              flag[index]  = add_z_vertex();
              set_z_vert(flag[index], _i, _j, _k);
            }
            //printf("3.3.2 %d\t", flag[index]);
          }
        }
      }
  index=0;
  for(tempi=0; tempi<_size_x*_size_y*_size_z; tempi++)    {
    if(flag[tempi]!=-1) {
      //printf("%d\n", flag[tempi]);
      index++;
    }
  }
  printf("sata:count = %d\n", index);
  free(flag);
}

bool test_face(int face)
{
  double A, B, C, D;
  switch(face)
    {
    case -1: case 1:
      A=_cube[0]; B = _cube[4] ;
      C = _cube[5] ;  D = _cube[1] ;
      break ;
    case -2: case 2:
      A = _cube[1] ;  B = _cube[5] ;
      C = _cube[6] ;  D = _cube[2] ;
      break ;
    case -3 : case 3 :
      A = _cube[2] ;  B = _cube[6] ;
      C = _cube[7] ;  D = _cube[3] ;
      break ;
    case -4 : case 4 :
      A = _cube[3] ;  B = _cube[7] ;
      C = _cube[4] ;  D = _cube[0] ;
      break ;
    case -5 : case 5 :
      A = _cube[0] ;  B = _cube[3] ;
      C = _cube[2] ;  D = _cube[1] ;
      break ;
    case -6 : case 6 :
      A = _cube[4] ;  B = _cube[7] ;
      C = _cube[6] ;  D = _cube[5] ;  break ;
    default :
      printf( "Invalid face code %d\n", face ) ;
      print_cube() ;  A = B = C = D = 0 ;
    };

  if( fabs( A*C - B*D ) < DBL_EPSILON )
    return (face >= 0 ? TRUE:FALSE );
  return (face * A * ( A*C - B*D ) >= 0 ? TRUE:FALSE) ;
}



bool test_interior(int s)   {
  double t, At=0, Bt=0, Ct=0, Dt=0, a, b;
  char  test =  0 ;
  char  edge = -1 ; // reference edge of the triangulation
  
  switch( _case )     {
  case  4 :
  case 10 :
    a = (_cube[4] - _cube[0]) * (_cube[6] - _cube[2]) -
      (_cube[7] - _cube[3] ) * ( _cube[5] - _cube[1]) ;
    b =  _cube[2] * ( _cube[4] - _cube[0] ) + _cube[0] * ( _cube[6] - _cube[2] )
      - _cube[1] * ( _cube[7] - _cube[3] ) - _cube[3] * ( _cube[5] - _cube[1] ) ;
    t = - b / (2*a) ;
    if( t<0 || t>1 ) return (s>0 ? TRUE:FALSE) ;
    
    At = _cube[0] + ( _cube[4] - _cube[0] ) * t ;
    Bt = _cube[3] + ( _cube[7] - _cube[3] ) * t ;
    Ct = _cube[2] + ( _cube[6] - _cube[2] ) * t ;
    Dt = _cube[1] + ( _cube[5] - _cube[1] ) * t ;
    break ;
  case  6 :
  case  7 :
  case 12 :
  case 13 :
    switch( _case )      {
    case  6 : edge = test6 [_config][2] ; break ;
    case  7 : edge = test7 [_config][4] ; break ;
    case 12 : edge = test12[_config][3] ; break ;
    case 13 : edge = tiling13_5_1[_config][_subconfig][0] ; break ;
    }
    switch( edge )  {
    case  0 :
      t  = _cube[0] / ( _cube[0] - _cube[1] ) ;
      At = 0 ;
      Bt = _cube[3] + ( _cube[2] - _cube[3] ) * t ;
      Ct = _cube[7] + ( _cube[6] - _cube[7] ) * t ;
      Dt = _cube[4] + ( _cube[5] - _cube[4] ) * t ;
      break ;
    case  1 :
      t  = _cube[1] / ( _cube[1] - _cube[2] ) ;
      At = 0 ;
      Bt = _cube[0] + ( _cube[3] - _cube[0] ) * t ;
      Ct = _cube[4] + ( _cube[7] - _cube[4] ) * t ;
      Dt = _cube[5] + ( _cube[6] - _cube[5] ) * t ;
      break ;
    case  2 :
      t  = _cube[2] / ( _cube[2] - _cube[3] ) ;
      At = 0 ;
      Bt = _cube[1] + ( _cube[0] - _cube[1] ) * t ;
      Ct = _cube[5] + ( _cube[4] - _cube[5] ) * t ;
      Dt = _cube[6] + ( _cube[7] - _cube[6] ) * t ;
      break ;
    case  3 :
      t  = _cube[3] / ( _cube[3] - _cube[0] ) ;
      At = 0 ;
      Bt = _cube[2] + ( _cube[1] - _cube[2] ) * t ;
      Ct = _cube[6] + ( _cube[5] - _cube[6] ) * t ;
      Dt = _cube[7] + ( _cube[4] - _cube[7] ) * t ;
      break ;
    case  4 :
      t  = _cube[4] / ( _cube[4] - _cube[5] ) ;
      At = 0 ;
      Bt = _cube[7] + ( _cube[6] - _cube[7] ) * t ;
      Ct = _cube[3] + ( _cube[2] - _cube[3] ) * t ;
      Dt = _cube[0] + ( _cube[1] - _cube[0] ) * t ;
      break ;
    case  5 :
      t  = _cube[5] / ( _cube[5] - _cube[6] ) ;
      At = 0 ;
      Bt = _cube[4] + ( _cube[7] - _cube[4] ) * t ;
      Ct = _cube[0] + ( _cube[3] - _cube[0] ) * t ;
      Dt = _cube[1] + ( _cube[2] - _cube[1] ) * t ;
      break ;
    case  6 :
      t  = _cube[6] / ( _cube[6] - _cube[7] ) ;
      At = 0 ;
      Bt = _cube[5] + ( _cube[4] - _cube[5] ) * t ;
      Ct = _cube[1] + ( _cube[0] - _cube[1] ) * t ;
      Dt = _cube[2] + ( _cube[3] - _cube[2] ) * t ;
      break ;
    case  7 :
      t  = _cube[7] / ( _cube[7] - _cube[4] ) ;
      At = 0 ;
      Bt = _cube[6] + ( _cube[5] - _cube[6] ) * t ;
      Ct = _cube[2] + ( _cube[1] - _cube[2] ) * t ;
      Dt = _cube[3] + ( _cube[0] - _cube[3] ) * t ;
      break ;
    case  8 :
      t  = _cube[0] / ( _cube[0] - _cube[4] ) ;
      At = 0 ;
      Bt = _cube[3] + ( _cube[7] - _cube[3] ) * t ;
      Ct = _cube[2] + ( _cube[6] - _cube[2] ) * t ;
      Dt = _cube[1] + ( _cube[5] - _cube[1] ) * t ;
      break ;
    case  9 :
      t  = _cube[1] / ( _cube[1] - _cube[5] ) ;
      At = 0 ;
      Bt = _cube[0] + ( _cube[4] - _cube[0] ) * t ;
      Ct = _cube[3] + ( _cube[7] - _cube[3] ) * t ;
      Dt = _cube[2] + ( _cube[6] - _cube[2] ) * t ;
      break ;
      
    case 10 :
      t  = _cube[2] / ( _cube[2] - _cube[6] ) ;
      At = 0 ;
      Bt = _cube[1] + ( _cube[5] - _cube[1] ) * t ;
      Ct = _cube[0] + ( _cube[4] - _cube[0] ) * t ;
      Dt = _cube[3] + ( _cube[7] - _cube[3] ) * t ;
      break ;
    case 11 :
      t  = _cube[3] / ( _cube[3] - _cube[7] ) ;
      At = 0 ;
      Bt = _cube[2] + ( _cube[6] - _cube[2] ) * t ;
      Ct = _cube[1] + ( _cube[5] - _cube[1] ) * t ;
      Dt = _cube[0] + ( _cube[4] - _cube[0] ) * t ;
      break ;
    default : printf( "Invalid edge %d\n", edge ) ;  print_cube() ;  break ;
    }
    break ;

  default : printf( "Invalid ambiguous case %d\n", _case ) ;  print_cube() ;  break ;
  }
  if( At >= 0 ) test ++ ;
  if( Bt >= 0 ) test += 2 ;
  if( Ct >= 0 ) test += 4 ;
  if( Dt >= 0 ) test += 8 ;
  switch( test )    {
  case  0 : return (s>0 ? TRUE:FALSE);
  case  1 : return (s>0 ? TRUE:FALSE);
  case  2 : return (s>0 ? TRUE:FALSE) ;
  case  3 : return (s>0 ? TRUE:FALSE) ;
  case  4 : return (s>0 ? TRUE:FALSE) ;
  case  5 : if( At * Ct - Bt * Dt <  DBL_EPSILON ) return (s>0 ? TRUE:FALSE) ; break ;
  case  6 : return (s>0 ? TRUE:FALSE) ;
  case  7 : return (s<0 ? TRUE:FALSE) ;
  case  8 : return (s>0 ? TRUE:FALSE) ;
  case  9 : return (s>0 ? TRUE:FALSE) ;
  case 10 : if( At * Ct - Bt * Dt >= DBL_EPSILON ) return (s>0 ? TRUE:FALSE) ; break ;
  case 11 : return (s<0 ? TRUE:FALSE) ;
  case 12 : return (s>0 ? TRUE:FALSE) ;
  case 13 : return (s<0 ? TRUE:FALSE) ;
  case 14 : return (s<0 ? TRUE:FALSE) ;
  case 15 : return (s<0 ? TRUE:FALSE) ;
  }
  
  return (s<0 ? TRUE:FALSE) ;
}

void process_cube()
{
  int v12=-1;

  //perform Classic Marching cube.
  if(_originalMC)   {
    int nt=0;
    while(casesClassic[_lut_entry][3*nt] != -1) nt++;
    add_triangle(casesClassic[_lut_entry], nt, -1);
    return;
  }
  //marching cube 33 implemented
  
  _case   = cases[_lut_entry][0] ;
  _config = cases[_lut_entry][1] ;
  _subconfig = 0 ;

  switch(_case)  {
  case 0:
    break;
  case  1 :
    add_triangle(tiling1[_config], 1, v12) ;
    break ;
  case  2 :
    add_triangle(tiling2[_config], 2, v12) ;
    break ;
  case  3 :
    if(test_face(test3[_config]) )
      add_triangle(tiling3_2[_config], 4, v12) ; // 3.2
    else
      add_triangle(tiling3_1[_config], 2, v12) ; // 3.1
    break ;
    
  case  4 :
    if(test_interior(test4[_config]) )
      add_triangle(tiling4_1[_config], 2, v12) ; // 4.1.1
    else
      add_triangle(tiling4_2[_config], 6, v12) ; // 4.1.2
    break ;
    
  case  5 :
    add_triangle(tiling5[_config], 3, v12) ;
    break ;
    
  case  6 :
    if(test_face(test6[_config][0]) )
      add_triangle(tiling6_2[_config], 5, v12) ; // 6.2
    else  {
      if(test_interior(test6[_config][1]) )
        add_triangle(tiling6_1_1[_config], 3, v12) ; // 6.1.1
      else
        add_triangle(tiling6_1_2[_config], 7, v12) ; // 6.1.2
    }
    break ;
    
  case  7 :
    if( test_face(test7[_config][0] ) ) _subconfig +=  1 ;
    if( test_face(test7[_config][1] ) ) _subconfig +=  2 ;
    if( test_face(test7[_config][2] ) ) _subconfig +=  4 ;
    switch( _subconfig )  {
    case 0 :
      add_triangle( tiling7_1[_config], 3, v12) ; 
      break ;
    case 1 :
      add_triangle(tiling7_2[_config][0], 5, v12) ;
      break ;
    case 2 :
      add_triangle(tiling7_2[_config][1], 5, v12) ;
      break ;
    case 3 :
      v12 = add_c_vertex() ;
      add_triangle(tiling7_3[_config][0], 9, v12) ; 
      break ;
    case 4 :
      add_triangle(tiling7_2[_config][2], 5, v12) ;
      break ;
    case 5 :
      v12 = add_c_vertex() ;
      add_triangle(tiling7_3[_config][1], 9, v12) ; 
      break ;
    case 6 :
      v12 = add_c_vertex() ;
      add_triangle(tiling7_3[_config][2], 9, v12) ;
      break ;
    case 7 :
      if(test_interior(test7[_config][3]))
        add_triangle(tiling7_4_2[_config], 9, v12) ;
      else
        add_triangle(tiling7_4_1[_config], 5, v12) ;
      break ;
    };
    break ;
    
  case  8 :
    add_triangle(tiling8[_config], 2, v12) ;
    break ;
    
  case  9 :
    add_triangle(tiling9[_config], 4, v12) ;
    break ;
    
  case 10 :
    if( test_face(test10[_config][0]))	{
      if(test_face(test10[_config][1]))
        add_triangle(tiling10_1_1_[_config], 4, v12) ; // 10.1.1
      else   {
        v12 = add_c_vertex() ;
        add_triangle(tiling10_2[_config], 8, v12) ; // 10.2
      }
    }
    else   {
      if(test_face( test10[_config][1]))  {
        v12 = add_c_vertex() ;
        add_triangle(tiling10_2_[_config], 8, v12) ; // 10.2
      }
      else  {
        if( test_interior(test10[_config][2]))
          add_triangle(tiling10_1_1[_config], 4, v12) ; // 10.1.1
        else
          add_triangle(tiling10_1_2[_config], 8, v12) ; // 10.1.2
      }
    }
    break ;
    
  case 11 :
    add_triangle(tiling11[_config], 4, v12) ;
    break ;
    
  case 12 :
    if(test_face( test12[_config][0]))	{
      if(test_face( test12[_config][1]))
        add_triangle(tiling12_1_1_[_config], 4, v12) ; // 12.1.1
      else   {
        v12 = add_c_vertex() ;
        add_triangle(tiling12_2[_config], 8, v12) ; // 12.2
      }
    }
    else   {
      if( test_face( test12[_config][1]) )   {
        v12 = add_c_vertex() ;
        add_triangle( tiling12_2_[_config], 8, v12) ; // 12.2
      }
      else   {
        if( test_interior( test12[_config][2]) )
          add_triangle( tiling12_1_1[_config], 4, v12) ; // 12.1.1
        else
          add_triangle( tiling12_1_2[_config], 8, v12) ; // 12.1.2
      }
    }
    break ;
    
  case 13 :
    if( test_face( test13[_config][0] ) ) _subconfig +=  1 ;
    if( test_face( test13[_config][1] ) ) _subconfig +=  2 ;
    if( test_face( test13[_config][2] ) ) _subconfig +=  4 ;
    if( test_face( test13[_config][3] ) ) _subconfig +=  8 ;
    if( test_face( test13[_config][4] ) ) _subconfig += 16 ;
    if( test_face( test13[_config][5] ) ) _subconfig += 32 ;
    switch( subconfig13[_subconfig] )	{
    case 0 :/* 13.1 */
      add_triangle( tiling13_1[_config], 4, v12) ; 
      break ;
      
    case 1 :/* 13.2 */
      add_triangle( tiling13_2[_config][0], 6, v12) ;
      break ;
    case 2 :/* 13.2 */
      add_triangle( tiling13_2[_config][1], 6, v12) ;
      break ;
    case 3 :/* 13.2 */
      add_triangle( tiling13_2[_config][2], 6, v12) ;
      break ;
    case 4 :/* 13.2 */
      add_triangle(tiling13_2[_config][3], 6, v12) ;
      break ;
    case 5 :/* 13.2 */
      add_triangle(tiling13_2[_config][4], 6, v12) ; 
      break ;
    case 6 :/* 13.2 */
      add_triangle( tiling13_2[_config][5], 6, v12) ;
      break ;
      
    case 7 :/* 13.3 */
      v12 = add_c_vertex() ;
      add_triangle(tiling13_3[_config][0], 10, v12 ) ; 
      break ;
    case 8 :/* 13.3 */
      v12 = add_c_vertex() ;
      add_triangle(tiling13_3[_config][1], 10, v12) ;
      break ;
    case 9 :/* 13.3 */
      v12 = add_c_vertex() ;
      add_triangle( tiling13_3[_config][2], 10, v12) ; 
      break ;
    case 10 :/* 13.3 */
      v12 = add_c_vertex() ;
      add_triangle(tiling13_3[_config][3], 10, v12) ;
      break ;
    case 11 :/* 13.3 */
      v12 = add_c_vertex() ;
      add_triangle( tiling13_3[_config][4], 10, v12) ;
      break ;
    case 12 :/* 13.3 */
      v12 = add_c_vertex() ;
      add_triangle(tiling13_3[_config][5], 10, v12) ;
      break ;
    case 13 :/* 13.3 */
      v12 = add_c_vertex() ;
      add_triangle(tiling13_3[_config][6], 10, v12) ;
      break ;
    case 14 :/* 13.3 */
      v12 = add_c_vertex() ;
      add_triangle(tiling13_3[_config][7], 10, v12) ;
      break ;
    case 15 :/* 13.3 */
      v12 = add_c_vertex() ;
      add_triangle( tiling13_3[_config][8], 10, v12) ;
      break ;
    case 16 :/* 13.3 */
      v12 = add_c_vertex() ;
      add_triangle(tiling13_3[_config][9], 10, v12) ;
      break ;
    case 17 :/* 13.3 */
      v12 = add_c_vertex() ;
      add_triangle( tiling13_3[_config][10], 10, v12) ; 
      break ;
    case 18 :/* 13.3 */
      v12 = add_c_vertex() ;
      add_triangle( tiling13_3[_config][11], 10, v12) ; 
      break ;

    case 19 :/* 13.4 */
      v12 = add_c_vertex() ;
      add_triangle( tiling13_4[_config][0], 12, v12) ; 
      break ;
    case 20 :/* 13.4 */
      v12 = add_c_vertex() ;
      add_triangle( tiling13_4[_config][1], 12, v12) ;
      break ;
    case 21 :/* 13.4 */
      v12 = add_c_vertex() ;
      add_triangle( tiling13_4[_config][2], 12, v12) ; 
      break ;
    case 22 :/* 13.4 */
      v12 = add_c_vertex() ;
      add_triangle( tiling13_4[_config][3], 12, v12) ; 
      break ;

    case 23 :/* 13.5 */
      _subconfig = 0 ;
      if( test_interior( test13[_config][6] ) )
        add_triangle(tiling13_5_1[_config][0], 6, v12) ;
      else
        add_triangle(tiling13_5_2[_config][0], 10, v12) ;
      break ;
    case 24 :/* 13.5 */
      _subconfig = 1 ;
      if( test_interior( test13[_config][6] ) )
        add_triangle(tiling13_5_1[_config][1], 6, v12) ;
      else
        add_triangle( tiling13_5_2[_config][1], 10, v12) ;
      break ;
    case 25 :/* 13.5 */
      _subconfig = 2 ;
      if( test_interior( test13[_config][6] ) )
        add_triangle( tiling13_5_1[_config][2], 6, v12) ;
      else
        add_triangle( tiling13_5_2[_config][2], 10, v12) ;
      break ;
    case 26 :/* 13.5 */
      _subconfig = 3 ;
      if( test_interior( test13[_config][6] ) )
        add_triangle( tiling13_5_1[_config][3], 6, v12) ;
      else
        add_triangle( tiling13_5_2[_config][3], 10, v12) ;
      break ;

    case 27 :/* 13.3 */
      v12 = add_c_vertex() ;
      add_triangle( tiling13_3_[_config][0], 10, v12) ; 
      break ;
    case 28 :/* 13.3 */
      v12 = add_c_vertex() ;
      add_triangle( tiling13_3_[_config][1], 10, v12) ; 
      break ;
    case 29 :/* 13.3 */
      v12 = add_c_vertex() ;
      add_triangle( tiling13_3_[_config][2], 10, v12) ;
      break ;
    case 30 :/* 13.3 */
      v12 = add_c_vertex() ;
      add_triangle( tiling13_3_[_config][3], 10, v12) ; 
      break ;
    case 31 :/* 13.3 */
      v12 = add_c_vertex() ;
      add_triangle( tiling13_3_[_config][4], 10, v12) ; 
      break ;
    case 32 :/* 13.3 */
      v12 = add_c_vertex() ;
      add_triangle( tiling13_3_[_config][5], 10, v12) ; 
      break ;
    case 33 :/* 13.3 */
      v12 = add_c_vertex() ;
      add_triangle( tiling13_3_[_config][6], 10, v12) ; 
      break ;
    case 34 :/* 13.3 */
      v12 = add_c_vertex() ;
      add_triangle( tiling13_3_[_config][7], 10, v12) ; 
      break ;
    case 35 :/* 13.3 */
      v12 = add_c_vertex() ;
      add_triangle( tiling13_3_[_config][8], 10, v12) ;
      break ;
    case 36 :/* 13.3 */
      v12 = add_c_vertex() ;
      add_triangle( tiling13_3_[_config][9], 10, v12) ;
      break ;
    case 37 :/* 13.3 */
      v12 = add_c_vertex() ;
      add_triangle( tiling13_3_[_config][10], 10, v12) ; 
      break ;
    case 38 :/* 13.3 */
      v12 = add_c_vertex() ;
      add_triangle( tiling13_3_[_config][11], 10, v12) ;
      break ;

    case 39 :/* 13.2 */
      add_triangle( tiling13_2_[_config][0], 6 , v12) ; break ;
    case 40 :/* 13.2 */
      add_triangle( tiling13_2_[_config][1], 6 , v12) ; break ;
    case 41 :/* 13.2 */
      add_triangle( tiling13_2_[_config][2], 6, v12 ) ; break ;
    case 42 :/* 13.2 */
      add_triangle( tiling13_2_[_config][3], 6, v12 ) ; break ;
    case 43 :/* 13.2 */
      add_triangle( tiling13_2_[_config][4], 6, v12 ) ; break ;
    case 44 :/* 13.2 */
      add_triangle( tiling13_2_[_config][5], 6, v12 ) ; break ;
    case 45 :/* 13.1 */
      add_triangle( tiling13_1_[_config], 4, v12) ; break ;

    default :
      printf("Marching Cubes: Impossible case 13?\n" ) ;  print_cube() ;
    }
    break ;

  case 14 :
    add_triangle( tiling14[_config], 4, v12) ;
    break ;
  }; 
}  //don't know, why the matching brackets are not shown

void add_triangle(int* trig, int n, int v12)   {
  int tv[3];   //locally triangle vertex
  int t;
  Triangle *T;
  
  
  for(t = 0; t < 3*n; t++)     {
    switch( trig[t] )      {
    case  0 : tv[ t % 3 ] = get_x_vert( _i , _j , _k ) ; break ;
    case  1 : tv[ t % 3 ] = get_y_vert(_i+1, _j , _k ) ; break ;
    case  2 : tv[ t % 3 ] = get_x_vert( _i ,_j+1, _k ) ; break ;
    case  3 : tv[ t % 3 ] = get_y_vert( _i , _j , _k ) ; break ;
    case  4 : tv[ t % 3 ] = get_x_vert( _i , _j ,_k+1) ; break ;
    case  5 : tv[ t % 3 ] = get_y_vert(_i+1, _j ,_k+1) ; break ;
    case  6 : tv[ t % 3 ] = get_x_vert( _i ,_j+1,_k+1) ; break ;
    case  7 : tv[ t % 3 ] = get_y_vert( _i , _j ,_k+1) ; break ;
    case  8 : tv[ t % 3 ] = get_z_vert( _i , _j , _k ) ; break ;
    case  9 : tv[ t % 3 ] = get_z_vert(_i+1, _j , _k ) ; break ;
    case 10 : tv[ t % 3 ] = get_z_vert(_i+1,_j+1, _k ) ; break ;
    case 11 : tv[ t % 3 ] = get_z_vert( _i ,_j+1, _k ) ; break ;
    case 12 : tv[ t % 3 ] = v12 ; break ;
    default : break ;
    }
    if( tv[t%3] == -1 )   {
      printf("Marching Cubes: invalid triangle %d %d %d %d %d \n", _ntrigs+1, trig[t], _i, _j, _k) ;
      print_cube() ;
    }
    if( t%3 == 2 )   {
      
      if(_Ntrigs <= _ntrigs)  {
        _Ntrigs *=2;
        _triangles = (Triangle *) realloc(_triangles, _Ntrigs*sizeof(Triangle));
      }
      T = _triangles + _ntrigs;
      T->v1 = tv[0];
      T->v2 = tv[1];
      T->v3 = tv[2];
      _ntrigs++;
    }
  }
}

double get_x_grad(int i, int j, int k)
{
  if(i > 0)
    {
      if (i < _size_x - 1)
        return (get_data(i+1, j, k) - get_data(i-1, j, k)) / 2. ;
      else
        return get_data(i, j, k) - get_data(i-1, j, k) ;
    }
  else
    return get_data(i+1, j, k) - get_data(i, j, k) ;
}

double get_y_grad(int i, int j, int k)
{
  if(j > 0)
    {
      if (j < _size_y - 1)
        return (get_data(i, j+1, k) - get_data(i, j-1, k)) / 2 ;
      else
        return get_data(i, j, k) - get_data(i, j-1, k) ;
    }
  else
    return get_data(i, j+1, k) - get_data(i, j, k) ;
}

double get_z_grad(int i, int j, int k)
{
  if(k > 0)
    {
      if (k < _size_z - 1)
        return (get_data(i, j, k+1) - get_data(i, j, k-1)) / 2. ;
      else
        return get_data(i, j, k) - get_data(i, j, k-1) ;
    }
  else
    return get_data(i, j, k+1) - get_data(i, j, k) ;
}

void test_vertex_addition()
{
  if( _nverts >= _Nverts )
    {
      _Nverts *= 2;
      _vertices = (Vertex *)realloc(_vertices, _Nverts*sizeof(Vertex));
    }
}

int add_x_vertex()
{
  double x1, x2, x3;
  double u = (_cube[0]) / (_cube[0] - _cube[1]) ;
  test_vertex_addition() ;

  /* if(fabs(u)<1.e-5)    u=0.; */
  /* if(fabs(u)>(1. - 1.e-5))    u=1.; */

  _vertices[_nverts].x      = (double)(_i+u)*LL;
  _vertices[_nverts].y      = (double)( _j )*LL;
  _vertices[_nverts].z      = (double)( _k )*LL;

  x1 = _vertices[_nverts].nx = (1-u)*get_x_grad(_i,_j,_k) + u*get_x_grad(_i+1,_j,_k);
  x2 = _vertices[_nverts].ny = (1-u)*get_y_grad(_i,_j,_k) + u*get_y_grad(_i+1,_j,_k);
  x3 = _vertices[_nverts].nz = (1-u)*get_z_grad(_i,_j,_k) + u*get_z_grad(_i+1,_j,_k);

  u=(double) sqrt(x1 * x1 + x2 * x2 + x3 * x3);

  if(u>0)
    {
      _vertices[_nverts].nx /= u ;
      _vertices[_nverts].ny /= u ;
      _vertices[_nverts].nz /= u ;
    }
  _nverts++;

  return _nverts-1;
}

int add_y_vertex()
{
  double x1, x2, x3;
  double u = (_cube[0]) / (_cube[0] - _cube[3]) ;
  test_vertex_addition() ;

  /* if(fabs(u)<EPSILON)    u=0.; */
  /* if(fabs(u)>(1. - EPSILON))    u=1.; */

  _vertices[_nverts].x      = (double)( _i  )*LL;
  _vertices[_nverts].y      = (double)( _j+u)*LL ;
  _vertices[_nverts].z      = (double)( _k  )*LL;

  x1 = _vertices[_nverts].nx = (1-u)*get_x_grad(_i,_j,_k) + u*get_x_grad(_i,_j+1,_k);
  x2 = _vertices[_nverts].ny = (1-u)*get_y_grad(_i,_j,_k) + u*get_y_grad(_i,_j+1,_k);
  x3 = _vertices[_nverts].nz = (1-u)*get_z_grad(_i,_j,_k) + u*get_z_grad(_i,_j+1,_k);

  u=(double) sqrt(x1 * x1 + x2 * x2 + x3 * x3);
  if(u>0)
    {
      _vertices[_nverts].nx /= u ;
      _vertices[_nverts].ny /= u ;
      _vertices[_nverts].nz /= u ;
    }
  _nverts++;

  return _nverts-1;
}

int add_z_vertex()
{
  double x1, x2, x3;
  double u = (_cube[0]) / (_cube[0] - _cube[4]) ;
  test_vertex_addition() ;

  /* if(fabs(u)<EPSILON)    u=0.; */
  /* if(fabs(u)>(1. - EPSILON))    u=1.; */

  _vertices[_nverts].x      = (double)( _i )*LL;
  _vertices[_nverts].y      = (double)( _j )*LL ;
  _vertices[_nverts].z      = (double)(_k+u)*LL ;

  x1 = _vertices[_nverts].nx = (1-u)*get_x_grad(_i,_j,_k) + u*get_x_grad(_i,_j,_k+1);
  x2 = _vertices[_nverts].ny = (1-u)*get_y_grad(_i,_j,_k) + u*get_y_grad(_i,_j,_k+1);
  x3 = _vertices[_nverts].nz = (1-u)*get_z_grad(_i,_j,_k) + u*get_z_grad(_i,_j,_k+1);

  u=(double) sqrt(x1 * x1 + x2 * x2 + x3 * x3);
  if(u>0)
    {
      _vertices[_nverts].nx /= u ;
      _vertices[_nverts].ny /= u ;
      _vertices[_nverts].nz /= u ;
    }
  _nverts++;

  return _nverts-1;
}

int add_c_vertex()
{
  double u = 0 ;
  int   vid ;
  test_vertex_addition() ;
  Vertex *vert = _vertices + _nverts;

  //_vertices[_nverts].x=_vertices[_nverts].y=_vertices[_nverts].z=0;
  //_vertices[_nverts].nx=_vertices[_nverts].ny=_vertices[_nverts].nz=0;
  vert->x = vert->y = vert->z =  vert->nx = vert->ny = vert->nz = 0 ;

  vid = get_x_vert( _i , _j , _k ) ;
  if(vid != -1)
    {
      ++u;
      vert -> x += _vertices[vid].x;      vert -> y += _vertices[vid].y;
      vert -> z += _vertices[vid].z;      vert -> nx += _vertices[vid].nx;
      vert -> ny += _vertices[vid].ny;    vert -> nz += _vertices[vid].nz;
    }

  vid = get_y_vert(_i+1, _j , _k ) ;
  if(vid != -1)
    {
      ++u;
      vert -> x += _vertices[vid].x;      vert -> y += _vertices[vid].y;
      vert -> z += _vertices[vid].z;      vert -> nx += _vertices[vid].nx;
      vert -> ny += _vertices[vid].ny;    vert -> nz += _vertices[vid].nz;
    }

  vid = get_x_vert( _i ,_j+1, _k ) ;
  if(vid != -1)
    {
      ++u;
      vert -> x += _vertices[vid].x;      vert -> y += _vertices[vid].y;
      vert -> z += _vertices[vid].z;      vert -> nx += _vertices[vid].nx;
      vert -> ny += _vertices[vid].ny;    vert -> nz += _vertices[vid].nz;
    }
  vid = get_y_vert( _i , _j , _k ) ;
  if(vid != -1)
    {
      ++u;
      vert -> x += _vertices[vid].x;      vert -> y += _vertices[vid].y;
      vert -> z += _vertices[vid].z;      vert -> nx += _vertices[vid].nx;
      vert -> ny += _vertices[vid].ny;    vert -> nz += _vertices[vid].nz;
    }

  vid = get_x_vert( _i , _j ,_k+1) ;
  if(vid != -1)
    {
      ++u;
      vert -> x += _vertices[vid].x;      vert -> y += _vertices[vid].y;
      vert -> z += _vertices[vid].z;      vert -> nx += _vertices[vid].nx;
      vert -> ny += _vertices[vid].ny;    vert -> nz += _vertices[vid].nz;
    }

  vid = get_y_vert(_i+1, _j ,_k+1) ;
  if(vid != -1)
    {
      ++u;
      vert -> x += _vertices[vid].x;      vert -> y += _vertices[vid].y;
      vert -> z += _vertices[vid].z;      vert -> nx += _vertices[vid].nx;
      vert -> ny += _vertices[vid].ny;    vert -> nz += _vertices[vid].nz;
    }

  vid = get_x_vert( _i ,_j+1,_k+1) ;
  if(vid != -1)
    {
      ++u;
      vert -> x += _vertices[vid].x;      vert -> y += _vertices[vid].y;
      vert -> z += _vertices[vid].z;      vert -> nx += _vertices[vid].nx;
      vert -> ny += _vertices[vid].ny;    vert -> nz += _vertices[vid].nz;
    }
  vid = get_y_vert( _i , _j ,_k+1) ;
  if(vid != -1)
    {
      ++u;
      vert -> x += _vertices[vid].x;      vert -> y += _vertices[vid].y;
      vert -> z += _vertices[vid].z;      vert -> nx += _vertices[vid].nx;
      vert -> ny += _vertices[vid].ny;    vert -> nz += _vertices[vid].nz;
    }

  vid = get_z_vert( _i , _j , _k ) ;
  if(vid != -1)
    {
      ++u;
      vert -> x += _vertices[vid].x;      vert -> y += _vertices[vid].y;
      vert -> z += _vertices[vid].z;      vert -> nx += _vertices[vid].nx;
      vert -> ny += _vertices[vid].ny;    vert -> nz += _vertices[vid].nz;
    }

  vid = get_z_vert(_i+1, _j , _k ) ;
  if(vid != -1)
    {
      ++u;
      vert -> x += _vertices[vid].x;      vert -> y += _vertices[vid].y;
      vert -> z += _vertices[vid].z;      vert -> nx += _vertices[vid].nx;
      vert -> ny += _vertices[vid].ny;    vert -> nz += _vertices[vid].nz;
    }

  vid = get_z_vert(_i+1,_j+1, _k ) ;
  if(vid != -1)
    {
      ++u;
      vert -> x += _vertices[vid].x;      vert -> y += _vertices[vid].y;
      vert -> z += _vertices[vid].z;      vert -> nx += _vertices[vid].nx;
      vert -> ny += _vertices[vid].ny;    vert -> nz += _vertices[vid].nz;
    }

  vid = get_z_vert( _i ,_j+1, _k ) ;
  if(vid != -1)
    {
      ++u;
      vert -> x += _vertices[vid].x;      vert -> y += _vertices[vid].y;
      vert -> z += _vertices[vid].z;      vert -> nx += _vertices[vid].nx;
      vert -> ny += _vertices[vid].ny;    vert -> nz += _vertices[vid].nz;
    }


  vert->x  /= u ;
  vert->y  /= u ;
  vert->z  /= u ;

  u = (double) sqrt( vert->nx * vert->nx + vert->ny * vert->ny +vert->nz * vert->nz ) ;
  if( u > 0 )
    {
      vert->nx /= u ;
      vert->ny /= u ;
      vert->nz /= u ;
    }
  _nverts++;
  return _nverts-1 ;
}



double get_data(int i, int j, int k)
{
  return(data[k + (j * _size_z) + (i * _size_z * _size_y)]);
}

int get_x_vert(int i, int j, int k )
{
  return _x_verts[ k + j*_size_z + i*_size_z*_size_y] ;
}

int   get_y_vert(  int i,  int j,  int k )
{
  return _y_verts[ k + j*_size_z + i*_size_z*_size_y] ;
}

int   get_z_vert(  int i,  int j,  int k )
{
  return _z_verts[ k + j*_size_z + i*_size_z*_size_y] ;
}

void  set_x_vert(int val, int i, int j,  int k )
{
  _x_verts[k + j*_size_z + i*_size_z*_size_y] = val;
}

void  set_y_vert(int val, int i, int j,  int k )
{
  _y_verts[k + j*_size_z + i*_size_z*_size_y] = val;
}

void  set_z_vert(int val, int i, int j,  int k )
{
  _z_verts[k + j*_size_z + i*_size_z*_size_y] = val;
}

void declare_TV_TNV()
{
      //int size_safe=NT+1;  
        TV = (double ***)malloc(NT*sizeof(double**));
        TNV= (Vector **)malloc(NT*sizeof(Vector*));
        int i; int j;
        for ( i = 0; i< NT; i++) // to avoid segmentation fault, <=NT2 is used. Actually <NT2 creates an array of size NT2. Still I am getting the segmentation fault in initialize(), in last of while. This problem does not appear if we use an array TV[NT2][3][3] instead!!
        {
                TV[i] = (double **) malloc(3*sizeof(double *));
                TNV[i] = (Vector *) malloc(3*sizeof(Vector));
                
                for ( j = 0; j < 3; j++) 
                {
                        TV[i][j] = (double *)malloc(3*sizeof(double));
                }

       }
}

int store_triangle(int NN)
{
	NT=0;
  int ii,jj=0;
  for(ii=0; ii<_ntrigs; ii++)
        {
          if(TriangleArea(_triangles[ii])>0.0) {NT++;}
        }
  
  declare_TV_TNV();      
  
  FILE *fp; FILE *fp1;
  char hname[100];
  
  //int count = 0;
  
  if(_ntrigs>1)
    {
    	printf("storing trianglas\n");
      count++;
     int zero_tri=0;

      if(NN==NN_check){  sprintf(hname, "%strianles_iso_%lf_NN%d.dat", fout, iso,NN);
      fp=fopen(hname, "w");
      }
      fp1=fopen("zero_tri.dat","w");
      for(ii=0; ii<_ntrigs; ii++)
        {
          if(TriangleArea(_triangles[ii])>0.0)
            //fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
            {
              //printf("%d %d  %d  %d\n", ii, _triangles[ii].v1, _triangles[ii].v2, _triangles[ii].v3);
               if(NN==NN_check){fprintf(fp, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
                      _vertices[_triangles[ii].v1].x, _vertices[_triangles[ii].v1].y,
                      _vertices[_triangles[ii].v1].z, _vertices[_triangles[ii].v2].x,
                      _vertices[_triangles[ii].v2].y, _vertices[_triangles[ii].v2].z,
                      _vertices[_triangles[ii].v3].x, _vertices[_triangles[ii].v3].y,
                      _vertices[_triangles[ii].v3].z);}
                      
                      TV[jj][0][0]=_vertices[_triangles[ii].v1].x; TV[jj][0][1]=_vertices[_triangles[ii].v1].y; TV[jj][0][2]=_vertices[_triangles[ii].v1].z;
                      TV[jj][1][0]=_vertices[_triangles[ii].v2].x; TV[jj][1][1]=_vertices[_triangles[ii].v2].y; TV[jj][1][2]=_vertices[_triangles[ii].v2].z;
                      TV[jj][2][0]=_vertices[_triangles[ii].v3].x; TV[jj][2][1]=_vertices[_triangles[ii].v3].y; TV[jj][2][2]=_vertices[_triangles[ii].v3].z;
                      jj++;
                      
            }
            else{zero_tri++; fprintf(fp1, " %lf\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",TriangleArea(_triangles[ii]),
                      _vertices[_triangles[ii].v1].x, _vertices[_triangles[ii].v1].y,
                      _vertices[_triangles[ii].v1].z, _vertices[_triangles[ii].v2].x,
                      _vertices[_triangles[ii].v2].y, _vertices[_triangles[ii].v2].z,
                      _vertices[_triangles[ii].v3].x, _vertices[_triangles[ii].v3].y,
                      _vertices[_triangles[ii].v3].z);}
        }
       if(NN==NN_check){fclose(fp);}
       fclose(fp1);
        printf("\n No of triangles: _trigs=%d, NT1= %d, zero_tri=%d\n", _ntrigs,NT,zero_tri);
      return NT;
    }
  else
    return -1;
}

double compute_vol()
{
  int ii;
  double  meanx, meany, meanz, volume, area;
  double x1, x2, x3, y1, y2, y3, z1, z2, z3;
  double areax, areay, areaz;

  volume=0.;
  for(ii=0; ii<_ntrigs; ii++)
    {
      meanx = _vertices[_triangles[ii].v1].x +
        _vertices[_triangles[ii].v2].x + _vertices[_triangles[ii].v3].x;
      meany = _vertices[_triangles[ii].v1].y +
        _vertices[_triangles[ii].v2].y + _vertices[_triangles[ii].v3].y;
      meanz = _vertices[_triangles[ii].v1].z +
        _vertices[_triangles[ii].v2].z + _vertices[_triangles[ii].v3].z;

      x1 = _vertices[_triangles[ii].v1].x;
      y1 = _vertices[_triangles[ii].v1].y;
      z1 = _vertices[_triangles[ii].v1].z;

      x2 = _vertices[_triangles[ii].v2].x - x1;
      y2 = _vertices[_triangles[ii].v2].y - y1;
      z2 = _vertices[_triangles[ii].v2].z - z1;

      x3 = _vertices[_triangles[ii].v3].x - x1;
      y3 = _vertices[_triangles[ii].v3].y - y1;
      z3 = _vertices[_triangles[ii].v3].z - z1;

      area = TriangleArea(_triangles[ii]);
      if(area>0)
        {
          areax = 0.5 * (y2 * z3 - y3 * z2) / area;
          areay = 0.5 * (x3 * z2 - x2 * z3) / area;
          areaz = 0.5 * (x2 * y3 - x3 * y2) / area;
        }
      else
        {
          areax = 0.5 * (y2 * z3 - y3 * z2);
          areay = 0.5 * (x3 * z2 - x2 * z3);
          areaz = 0.5 * (x2 * y3 - x3 * y2);
        }
      area *= (-1.0/9.0)*((meanx * areax) + (meany * areay)
                          + (meanz * areaz));
      volume += area;

    }
  return((volume));
}

double TriangleArea(Triangle T)
{
  double x1, y1, z1, x2, y2, z2, x3, y3, z3;
  double area;

  x1 = (double)_vertices[T.v1].x;
  y1 = (double)_vertices[T.v1].y;
  z1 = (double)_vertices[T.v1].z;

  x2 = (double)_vertices[T.v2].x - x1;
  y2 = (double)_vertices[T.v2].y - y1;
  z2 = (double)_vertices[T.v2].z - z1;

  x3 = (double)_vertices[T.v3].x - x1;
  y3 = (double)_vertices[T.v3].y - y1;
  z3 = (double)_vertices[T.v3].z - z1;

  //printf("%e\n", area);
  area = (double)0.5 * sqrt(pow(y2 * z3 - y3 * z2, (double)2.) +
                            pow(x3 * z2 - x2 * z3, (double)2.) +
                            pow(x2 * y3 - x3 * y2, (double)2.));

  return ((double)area);
}

double compute_area()
{
  int ii;
  double total_area=0.;

  for(ii=0; ii<_ntrigs; ii++)
    {
      total_area += TriangleArea(_triangles[ii]);
    }
  return(total_area);
}

double common_edge(Triangle A, Triangle B)
{
  int ii, jj, kk=0;
  double length;
  double MX[3], MY[3], MZ[3];
  Vertex VA[3], VB[3];

  VA[0] = _vertices[A.v1];
  VA[1] = _vertices[A.v2];
  VA[2] = _vertices[A.v3];

  VB[0] = _vertices[B.v1];
  VB[1] = _vertices[B.v2];
  VB[2] = _vertices[B.v3];

  for(ii=0; ii<3; ii++)
    for(jj=0; jj<3; jj++)
      {
        if((fabs(VA[ii].x - VB[jj].x) < DBL_EPSILON) &&
           (fabs(VA[ii].y - VB[jj].y) < DBL_EPSILON) &&
           (fabs(VA[ii].z - VB[jj].z) < DBL_EPSILON))
          {
            MX[kk]=(double) VA[ii].x;
            MY[kk]=(double) VA[ii].y;
            MZ[kk]=(double) VA[ii].z;
            kk++;
          }
      }

  if(kk==2)
    length = sqrt(pow((MX[0] - MX[1]), 2.) +
                  pow((MY[0] - MY[1]), 2.) +
                  pow((MZ[0] - MZ[1]), 2.));
  else length = 0;
  return length;
}


/* ========================================================= */

int trind;
int make_hash(int *hashval, htree *A, int i, int flag) //htree defined in MarchingCube.h
{
  htree *tangle;

  if(i == (2*MAXLENGTH-1))
    {
      if(hashval[i]==0)
        {
          if(A->next0 == NULL)
            {
              tangle = (htree *)malloc(sizeof(htree));
              tangle->tr = trind;
              A->next0 = tangle;
              return 0;
            }
          else
            {
              tangle = A->next0;
              trind = tangle->tr;
              free(tangle);
              A->next0 = NULL;
              return 1;
            }

        }
      else
        {

          if(A->next1 == NULL)
            {
              tangle = (htree *)malloc(sizeof(htree));
              tangle->tr = trind;
              A->next1 = tangle;
              return 0;
            }
          else
            {
              tangle = A->next1;
              trind = tangle->tr;
              free(tangle);
              A->next1 = NULL;
              return 1;
            }
        }
    }

  else
    {
      if(hashval[i]==0)
        {
          if(A->next0 == NULL)
            {
              htree *B;
              B = (htree *)malloc(sizeof(htree));
              B->next0=NULL;
              B->next1=NULL;
              B->tr = A->tr;
              A->next0 = B;
              return make_hash(hashval,B,++i,flag);
            }
          else
            {
              flag = make_hash(hashval,A->next0,++i,flag);
              if(flag == 1)
                {
                  A->next0 = NULL;
                  if(A->next1 ==NULL)
                    {
                      free(A);
                      return 1;
                    }
                  else
                    return 0;
                }
              else if(flag == 0)
                return flag;
            }
        }
      else
        {
          if(A->next1 == NULL)
            {
              htree *B;
              B = (htree *)malloc(sizeof(htree));
              B->next0=NULL;
              B->next1=NULL;
              B->tr = A->tr;
              A->next1 = B;
              return make_hash(hashval,B,++i,flag);
            }
          else
            {
              flag = make_hash(hashval,A->next1,++i,flag);
              if(flag == 1)
                {
                  A->next1 = NULL;
                  if(A->next0 ==NULL)
                    {
                      free(A);
                      return 1;
                    }
                  else
                    return 0;
                }

              else if(flag == 0)
                return flag;
            }
        }
    }
  return -1;
}

void hash(unsigned int v1, unsigned int v2, int b1[])
{
  int i=0, j;
  for(j=0;j<2*MAXLENGTH;j++)
    b1[j]= 0;
  do
    {
      b1[i++]=v1%2;
      v1 = v1/2;
    }
  while(v1 > 0);

  i = MAXLENGTH;
  do
    {
      b1[i++]=v2%2;
      v2 = v2/2;
    }
  while(v2 > 0);
}

void printB(int b[])
{
  int i;
  for(i=0; i<2*MAXLENGTH; i++)
    printf("%d", b[i]);
  printf("\n");
}

//-------curvature computation based on Trimesh-2 ----------

Vector CROSSPRODUCT(Vector V1, Vector V2)
{
  Vector V;
  V.x = (V1.y) * (V2.z) - (V1.z) * (V2.y);
  V.y = (V1.z) * (V2.x) - (V1.x) * (V2.z);
  V.z = (V1.x) * (V2.y) - (V1.y) * (V2.x);

  return V;
}

double DOTPRODUCT(Vector V1, Vector V2)
{
  double dpdt;
  dpdt = (V1.x)*(V2.x) + (V1.y)*(V2.y) + (V1.z)*(V2.z);
  return dpdt;
}

//Rotate a coordinate system to be perpendicular to the given normal

void rot_coord_sys(Vector old_u, Vector old_v, Vector new_norm,
                   Vector *new_u1, Vector *new_v1)
{
  double ndot;
  Vector old_norm, perp_old, dperp;
  Vector new_u, new_v;

  new_u = old_u;
  new_v = old_v;
  old_norm = CROSSPRODUCT(old_u, old_v);
  ndot = DOTPRODUCT(old_norm, new_norm);

  if(ndot<=-1)
    {
      new_u.x=-new_u.x;  new_u.y=-new_u.y;  new_u.z=-new_u.z;
      new_v.x=-new_v.x;  new_v.y=-new_v.y;  new_v.z=-new_v.z;

      *new_u1 = new_u;
      *new_v1 = new_v;
      return;
    }

  perp_old.x = new_norm.x - ndot * old_norm.x;
  perp_old.y = new_norm.y - ndot * old_norm.y;
  perp_old.z = new_norm.z - ndot * old_norm.z;

  dperp.x = 1. / (1+ndot) * (old_norm.x + new_norm.x);
  dperp.y = 1. / (1+ndot) * (old_norm.y + new_norm.y);
  dperp.z = 1. / (1+ndot) * (old_norm.z + new_norm.z);

  new_u.x -= dperp.x * DOTPRODUCT(new_u, perp_old);
  new_u.y -= dperp.y * DOTPRODUCT(new_u, perp_old);
  new_u.z -= dperp.z * DOTPRODUCT(new_u, perp_old);

  new_v.x -= dperp.x * DOTPRODUCT(new_v, perp_old);
  new_v.y -= dperp.y * DOTPRODUCT(new_v, perp_old);
  new_v.z -= dperp.z * DOTPRODUCT(new_v, perp_old);


  *new_u1 = new_u;
  *new_v1 = new_v;
}

void proj_curv(Vector old_u , Vector old_v, double old_ku,
               double old_kuv, double old_kv, Vector new_u,
               Vector new_v, double *new_ku1, double *new_kuv1,
               double *new_kv1)
{
  Vector r_new_u, r_new_v;
  double new_ku, new_kuv, new_kv;
  double u1, u2, v1, v2;

  rot_coord_sys(new_u, new_v, CROSSPRODUCT(old_u, old_v), &r_new_u, &r_new_v);

  u1 = DOTPRODUCT(r_new_u, old_u);
  v1 = DOTPRODUCT(r_new_u, old_v);
  u2 = DOTPRODUCT(r_new_v, old_u);
  v2 = DOTPRODUCT(r_new_v, old_v);

  new_ku = (old_ku * u1*u1) + old_kuv * (2.0  * u1*v1) + old_kv * v1*v1;
  new_ku = (fabs(new_ku)<DBL_EPSILON)? 0:new_ku;
  new_kuv = old_ku * u1*u2 + old_kuv * (u1*v2 + u2*v1) + old_kv * v1*v2;
  new_kuv = (fabs(new_kuv)<DBL_EPSILON)? 0:new_kuv;
  new_kv  = old_ku * u2*u2 + old_kuv * (2.0  * u2*v2) + old_kv * v2*v2;
  new_kv = (fabs(new_kv)<DBL_EPSILON)? 0:new_kv;

  *new_ku1 = new_ku;
  *new_kuv1 = new_kuv;
  *new_kv1 = new_kv;
}

Vector Normalise(Vector V)
{
  double norm;
  norm = sqrt(V.x * V.x + V.y * V.y + V.z * V.z);

  V.x /= norm;
  V.y /= norm;
  V.z /= norm;

  return V;
}

double len2(Vector a)
{
  return (a.x * a.x + a.y * a.y + a.z * a.z);
}
double len(Vector a)
{
  return sqrt(len2(a));
}

void VectorNormal()
{
  int ii;

  Normals = (Vector *)calloc(_nverts, sizeof(Vector));
  for(ii = 0; ii<_nverts; ii++)
    {
      Normals[ii].x = _vertices[ii].nx;
      Normals[ii].y = _vertices[ii].ny;
      Normals[ii].z = _vertices[ii].nz;

      Normals[ii] = Normalise(Normals[ii]);
    }
}

double compute_curvature()
{
  int ii, jj, vj, VV[3];
  Vector  e[3], t, n, b;
  double m[3], w[3][3], u, v,dnu, dnv, diag[3], c1, c2, wt, c12;
  Vector cornerareas[_ntrigs], dn;
  double curv12[_nverts], area, l2[3], ew[3], ewscale, pointareas[_nverts];
  double curv1[_nverts], curv2[_nverts], sum=0;

  VectorNormal();


  for(ii=0; ii<_nverts; ii++)
    {
      pointareas[ii]=0;
      curv1[ii]=0;curv12[ii]=0;curv2[ii]=0;
    }

  for(ii=0; ii<_ntrigs; ii++)
    {
      cornerareas[ii].x = 0;
      cornerareas[ii].y = 0;
      cornerareas[ii].z = 0;
    }
  //compute corner areas

  for(ii = 0; ii < _ntrigs; ii++)
    {
      if(TriangleArea(_triangles[ii])) {
        e[0].x = _vertices[_triangles[ii].v3].x - _vertices[_triangles[ii].v2].x;
        e[0].y = _vertices[_triangles[ii].v3].y - _vertices[_triangles[ii].v2].y;
        e[0].z = _vertices[_triangles[ii].v3].z - _vertices[_triangles[ii].v2].z;

        e[1].x = _vertices[_triangles[ii].v1].x - _vertices[_triangles[ii].v3].x;
        e[1].y = _vertices[_triangles[ii].v1].y - _vertices[_triangles[ii].v3].y;
        e[1].z = _vertices[_triangles[ii].v1].z - _vertices[_triangles[ii].v3].z;

        e[2].x = _vertices[_triangles[ii].v2].x - _vertices[_triangles[ii].v1].x;
        e[2].y = _vertices[_triangles[ii].v2].y - _vertices[_triangles[ii].v1].y;
        e[2].z = _vertices[_triangles[ii].v2].z - _vertices[_triangles[ii].v1].z;

        area = 0.5 * len(CROSSPRODUCT(e[0], e[1]));
        //printf("%e\t%e\n", area, TriangleArea(_triangles[ii]));
        l2[0] = len2(e[0]);      l2[1] = len2(e[1]);      l2[2] = len2(e[2]);
        //printf("%e\t%e\t%e\n", l2[0], l2[1], l2[2]);
        ew[0] = l2[0] * (l2[1] + l2[2] - l2[0]);
        ew[1] = l2[1] * (l2[2] + l2[0] - l2[1]);
        ew[2] = l2[2] * (l2[0] + l2[1] - l2[2]);

        if(ew[0]<= 0.)
          {
            cornerareas[ii].y = -0.25 * l2[2] * area / DOTPRODUCT(e[0], e[2]);
            cornerareas[ii].z = -0.25 * l2[1] * area / DOTPRODUCT(e[0], e[1]);
            cornerareas[ii].x = area - cornerareas[ii].z - cornerareas[ii].y;
          }
        else if(ew[1] <= 0.)
          {
            cornerareas[ii].z = -0.25 * l2[0] * area / DOTPRODUCT(e[0], e[1]);
            cornerareas[ii].x = -0.25 * l2[2] * area / DOTPRODUCT(e[2], e[1]);
            cornerareas[ii].y = area - cornerareas[ii].x - cornerareas[ii].z;
          }
        else if(ew[2] <= 0.)
          {
            cornerareas[ii].x = -0.25 * l2[1] * area / DOTPRODUCT(e[2], e[1]);
            cornerareas[ii].y = -0.25 * l2[0] * area / DOTPRODUCT(e[2], e[0]);
            cornerareas[ii].z = area - cornerareas[ii].x - cornerareas[ii].y;
          }
        else
          {
            ewscale = 0.5 * area / (ew[0] + ew[1] + ew[2]);
            cornerareas[ii].x = ewscale * (ew[(0+1)%3] + ew[(0+2)%3]);
            cornerareas[ii].y = ewscale * (ew[(1+1)%3] + ew[(1+2)%3]);
            cornerareas[ii].z = ewscale * (ew[(2+1)%3] + ew[(2+2)%3]);
          }
        pointareas[_triangles[ii].v1] += cornerareas[ii].x;
        pointareas[_triangles[ii].v2] += cornerareas[ii].y;
        pointareas[_triangles[ii].v3] += cornerareas[ii].z;
      }
      else  {
        cornerareas[ii].x = 0;
        cornerareas[ii].y = 0;
        cornerareas[ii].z = 0;
        pointareas[_triangles[ii].v1] += 0;
        pointareas[_triangles[ii].v2] += 0;
        pointareas[_triangles[ii].v3] += 0;
      }
    }

  pdir1 = (Vector *)calloc(_nverts, sizeof(Vector)); /* direction of the first principal curvature */
  pdir2 = (Vector *)calloc(_nverts, sizeof(Vector)); /* direction of the second principal curvature */

  for(ii=0; ii<_ntrigs; ii++)    {
    if(TriangleArea(_triangles[ii])>DBL_EPSILON)    {
      pdir1[_triangles[ii].v1].x = _vertices[_triangles[ii].v2].x - _vertices[_triangles[ii].v1].x;
      pdir1[_triangles[ii].v1].y = _vertices[_triangles[ii].v2].y - _vertices[_triangles[ii].v1].y;
      pdir1[_triangles[ii].v1].z = _vertices[_triangles[ii].v2].z - _vertices[_triangles[ii].v1].z;

      pdir1[_triangles[ii].v2].x = _vertices[_triangles[ii].v3].x - _vertices[_triangles[ii].v2].x;
      pdir1[_triangles[ii].v2].y = _vertices[_triangles[ii].v3].y - _vertices[_triangles[ii].v2].y;
      pdir1[_triangles[ii].v2].z = _vertices[_triangles[ii].v3].z - _vertices[_triangles[ii].v2].z;

      pdir1[_triangles[ii].v3].x = _vertices[_triangles[ii].v1].x - _vertices[_triangles[ii].v3].x;
      pdir1[_triangles[ii].v3].y = _vertices[_triangles[ii].v1].y - _vertices[_triangles[ii].v3].y;
      pdir1[_triangles[ii].v3].z = _vertices[_triangles[ii].v1].z - _vertices[_triangles[ii].v3].z;
    }
  }

  for (ii = 0; ii < _nverts; ii++)    {
    pdir1[ii] = CROSSPRODUCT(pdir1[ii], Normals[ii]);
    pdir1[ii] = Normalise(pdir1[ii]);
    pdir2[ii] = CROSSPRODUCT(Normals[ii], pdir1[ii]);
    pdir2[ii] = Normalise(pdir2[ii]);
  }

  // Compute curvature per-face
  for (ii = 0; ii < _ntrigs; ii++)    {
    if(TriangleArea(_triangles[ii])>DBL_EPSILON) {
      e[0].x = _vertices[_triangles[ii].v3].x - _vertices[_triangles[ii].v2].x;
      e[0].y = _vertices[_triangles[ii].v3].y - _vertices[_triangles[ii].v2].y;
      e[0].z = _vertices[_triangles[ii].v3].z - _vertices[_triangles[ii].v2].z;

      e[1].x = _vertices[_triangles[ii].v1].x - _vertices[_triangles[ii].v3].x;
      e[1].y = _vertices[_triangles[ii].v1].y - _vertices[_triangles[ii].v3].y;
      e[1].z = _vertices[_triangles[ii].v1].z - _vertices[_triangles[ii].v3].z;

      e[2].x = _vertices[_triangles[ii].v2].x - _vertices[_triangles[ii].v1].x;
      e[2].y = _vertices[_triangles[ii].v2].y - _vertices[_triangles[ii].v1].y;
      e[2].z = _vertices[_triangles[ii].v2].z - _vertices[_triangles[ii].v1].z;

      t = e[0];
      t = Normalise(t);
      n=CROSSPRODUCT(e[0], e[1]);
      b=CROSSPRODUCT(n, t);
      b = Normalise(b);

      m[0] = 0;         m[1] = 0;         m[2] = 0;
      w[0][0] = 0;      w[1][0] = 0;      w[2][0] = 0;
      w[0][1] = 0;      w[1][1] = 0;      w[2][1] = 0;
      w[0][2] = 0;      w[1][2] = 0;      w[2][2] = 0;

      VV[0]=_triangles[ii].v1;
      VV[1]=_triangles[ii].v2;
      VV[2]=_triangles[ii].v3;

      for(jj = 0; jj < 3; jj++)    {
        u = DOTPRODUCT(e[jj], t);
        v = DOTPRODUCT(e[jj], b);
        w[0][0] += u*u;
        w[0][1] += u*v;
        w[2][2] += v*v;

        dn.x = Normals[VV[PREV(jj)]].x - Normals[VV[NEXT(jj)]].x;
        dn.y = Normals[VV[PREV(jj)]].y - Normals[VV[NEXT(jj)]].y;
        dn.z = Normals[VV[PREV(jj)]].z - Normals[VV[NEXT(jj)]].z;

        dnu = DOTPRODUCT(dn, t);
        dnv = DOTPRODUCT(dn, b);
        //printf("dnu %e\t%e\n", dnu, dnv);    
        m[0] += dnu*u;  m[1] += dnu*v + dnv*u;  m[2] += dnv*v;

        //printf("m[%d] %e\t%e\t%e\n", jj, m[0], m[1], m[2]);
      }
      w[1][1] = w[0][0] + w[2][2];
      w[1][2] = w[0][1];

      diag[0]=0;diag[1]=0;diag[2]=0;
      if(!ldltdc(w, diag))
        continue;
      ldltsl(w, diag, m, m);

      VV[0]=_triangles[ii].v1;
      VV[1]=_triangles[ii].v2;
      VV[2]=_triangles[ii].v3;

      vj=VV[0];
      proj_curv(t, b, m[0], m[1], m[2], pdir1[vj], pdir2[vj], &c1, &c12, &c2);

      wt = cornerareas[ii].x / pointareas[vj];
      curv1[vj]  += wt * c1;
      curv12[vj] += wt * c12;
      curv2[vj]  += wt * c2;

      vj=VV[1];
      proj_curv(t, b, m[0], m[1], m[2], pdir1[vj], pdir2[vj], &c1, &c12, &c2);
      wt = cornerareas[ii].y / pointareas[vj];
      curv1[vj]  += wt * c1;
      curv12[vj] += wt * c12;
      curv2[vj]  += wt * c2;

      vj=VV[2];
      proj_curv(t, b, m[0], m[1], m[2], pdir1[vj], pdir2[vj], &c1, &c12, &c2);
      wt = cornerareas[ii].z / pointareas[vj];
      curv1[vj]  += wt * c1;
      curv12[vj] += wt * c12;
      curv2[vj]  += wt * c2;
    }
  }
  for (ii = 0; ii<_nverts; ii++)  {
    diagonalize_curv(pdir1[ii], pdir2[ii], curv1[ii], curv12[ii], curv2[ii], Normals[ii], &pdir1[ii], &pdir2[ii], &curv1[ii], &curv2[ii]);
    sum += 0.5*(curv1[ii] + curv2[ii])*pointareas[ii];
  }

  free(Normals);

  free(pdir1);
  free(pdir2);

  return fabs(sum);
}

bool ldltdc(double A[3][3], double rdiag[3])
{
  double v[3], sum;
  int i, j, k;

  for (i =0; i< 3; i++)    {
    for (k = 0; k < i; k++)
      v[k] = A[i][k] * rdiag[k];
    for (j = i; j < 3; j++)   {
      sum =A[i][j];
      for (k = 0; k < i; k++)
        sum -= v[k] * A[j][k];
      if(i == j)  {
        if(sum<=0)
          return FALSE;
        rdiag[i]=1./sum;
      }
      else  {
        A[j][i]=sum;
      }

    }
  }
  return TRUE;
}

void ldltsl(double A[3][3], double rdiag[3], double B[3], double x[3])
{
  int i, k;
  double sum;

  for (i = 0; i < 3; i++)  {
    sum =  B[i];
    for (k = 0; k < i; k++)
      sum -= A[i][k] * x[k];
    x[i] = sum * rdiag[i];
  }
  for (i = 2; i >= 0; i--)  {
    sum =  0;
    for (k = i + 1; k < 3; k++)
      sum += A[k][i] * x[k];
    x[i] -= sum * rdiag[i];
  }
  for(i=0; i<3; i++)
    x[i] = (fabs(x[i])<DBL_EPSILON)? 0:x[i];
  //printf("x = %e\t%e\t%e\n", x[0], x[1], x[2]);
}

void diagonalize_curv(Vector old_u, Vector old_v, double ku, double kuv, double kv, Vector new_norm,
                      Vector *pdir11, Vector *pdir21, double *k11, double *k21)
{

  Vector r_old_u, r_old_v;
  double h, c = 1, s = 0, tt = 0;
  Vector pdir1, pdir2;
  double k1, k2;

  rot_coord_sys(old_u, old_v, new_norm, &r_old_u, &r_old_v);

  pdir1=*pdir11;
  pdir2=*pdir21;
  //printf(" .... %e\t%e\t%e\n", ku, kuv, kv);
  if (kuv!=0.)   {
    //Jacobi rotation to diagonalize
    h = 0.5 * (kv - ku)/kuv;
    tt = (h < 0.)? 1./(h - sqrt(1 + h*h)): 1./(h + sqrt(1 + h*h));
    //tt = (h + sqrt(1 + h*h));
    c = 1./sqrt(1 + tt*tt);
    s = tt * c;
  }
  k1 = ku - tt * kuv;
  k2 = kv + tt * kuv;

  if(fabs(k1) >= fabs(k2))  {
    pdir1.x = c*r_old_u.x - s*r_old_v.x;
    pdir1.y = c*r_old_u.y - s*r_old_v.y;
    pdir1.z = c*r_old_u.z - s*r_old_v.z;
  }
  else   {
    swap(&k1, &k2);
    pdir1.x = s*r_old_u.x + c*r_old_v.x;
    pdir1.y = s*r_old_u.y + c*r_old_v.y;
    pdir1.z = s*r_old_u.z + c*r_old_v.z;
  }

  *pdir11 = pdir1;
  *pdir21 = pdir2;
  *k11 = k1; *k21 = k2;
}

void swap(double *k1, double *k2)  {
  double temp;
  temp = *k1;
  *k1 = *k2;
  *k2 = temp;
}

