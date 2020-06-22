#ifndef MARCHINGCUBE_H_
#define MARCHINGCUBE_H_
# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <unistd.h>
# include <float.h>
# include <assert.h>
# include <string.h>
# include <malloc.h>

typedef enum{FALSE,TRUE}bool;

typedef struct {
  double  x,  y,  z ;  /**< Vertex coordinates */
  double nx, ny, nz ;  /**< Vertex normal */
} Vertex ;

typedef struct {
  int v1,v2,v3 ;  /**< Triangle vertices */
} Triangle ;

typedef struct node {
  int tr;
  struct node *next0,*next1;
} htree;

typedef struct
{
  double  x,  y,  z ;  /**< Vertex coordinates */
} Vector ;


int *_clusters;
int _size_x, _size_y, _size_z;
double *data;
double iso, LL;
char fout[100]; //to save triangles for a cluster
char fname[100]; //to save all the shapefinders for a rho_th:to distinguish between inputfile
char fshape[100]; //to save all the shapefinders for a rho_th
char fcluster_stat[100]; //optional, to save the cluster statistics at that threshold

double *** TV; Vector ** TNV;
int NT;
double vol, area, imc, genus;
double imc1, imc2;
int no_of_tri;
//int NT1=0; //for genus calculation; NE and NV are in cal_genus.c //NT1 is calculated in store_triangle() in Merchincube.c
//we name it NT1 because NT is also used globaly in some functions

int NN_check; //check points and triangles of this cluster. set in main.c.

void run();
int store_triangle();
double compute_vol();
double compute_area();
double compute_curvature();
double compute_genus_imc();
int Drive_cluster(char *);

void print_cube();
void print_vertex(int);
void init_temps();
void init_all();
void clean_temps();
void clean_all();
void compute_intersection_points();
bool test_face(int);

void process_cube();
void add_triangle(int*, int, int);
int add_y_vertex();
double get_data(int, int, int);
void  set_x_vert(int, int, int, int);
void  set_y_vert(int, int, int, int);
void  set_z_vert(int, int, int, int);
void swap(double*, double*);
double TriangleArea(Triangle);
Vector Normalise(Vector);


Vector CROSSPRODUCT(Vector, Vector);
double DOTPRODUCT(Vector, Vector);
void VectorNormal();
bool ldltdc(double A[3][3], double rdiag[3]);
void ldltsl(double A[3][3], double rdiag[3] , double B[3], double x[3]);
void diagonalize_curv(Vector, Vector, double, double, double, Vector, Vector*, Vector*, double*, double*);
int add_y_vertex();
int add_x_vertex();
int add_z_vertex();

int interior_ambiguity(int, int);
int interior_ambiguity_verification(int);
int add_c_vertex();

int get_x_vert(int, int, int);
int get_y_vert(int, int, int);
int get_z_vert(int, int, int);



double len2(Vector );
double len(Vector);
double compute_curvature();
void printsize();

int bad_NE;
int count_th;




#endif /* MARCHINGCUBE_H_ */
