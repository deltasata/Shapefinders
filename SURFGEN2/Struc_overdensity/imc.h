#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

double NE1;
double NE;



  enum nodeColor {
        RED,
        BLACK
  };

/*  
typedef struct
{
  double  x,  y,  z ;  // Vertex coordinates 
} Vector ;

typedef struct {
  double  x,  y,  z ;  //< Vertex coordinates 
  double nx, ny, nz ;  //< Vertex normal 
} Vertex ;

double *** TV;
Vector ** TNV;*/

  struct VertexNode {
        int color;
        //int conn_vert_no;
        int conn_tri_no;
        
       struct ConnectedVertNode *root_cvert;
        //int ConnectedVert[100]; //conn_vert_no_max=100 is a prior constant: Ver_no=3*tn+vn; tn=Ver_no/3; vn=Ver_no%3;
        double x,y,z;
        Vector nv;
        int shared_tri[100];  // Ver_no=3*tn+vn; tn=Ver_no/3; vn=Ver_no%3;
        //struct SharedTriNode *root_stri;
        struct VertexNode *link[2];
  };

  struct VertexNode *root; //root of the tree which stores Vertices
  
 
 
   struct ConnectedVertNode {
        int color; 
        int conn_vert_no; // Ver_no=3*tn+vn; tn=Ver_no/3; vn=Ver_no%3;
        struct ConnectedVertNode *link[2];
  };
  
  /*struct SharedTriNode {
        int color; 
        int shared_tri_no; // Ver_no=3*tn+vn; tn=Ver_no/3; vn=Ver_no%3;
        struct SharedTriNode *link[2];
  };*/
  
  
  
  void fixup_VertexNode(struct VertexNode **, int *, int );
  void insert_conn_vert(struct VertexNode * ,int);
  double get_ep_theta(int, int);

 void free_VertexNode(struct VertexNode *);
