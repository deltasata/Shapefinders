# include "MarchingCube.h"
#include "imc.h"
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>


void fixup_VertexNode(struct VertexNode *stack[98], int dir[98], int ht)
{

	struct VertexNode *xPtr, *yPtr;
        while ((ht >= 3) && (stack[ht - 1]->color == RED)) {
                if (dir[ht - 2] == 0) {
                        yPtr = stack[ht - 2]->link[1];
                        if (yPtr != NULL && yPtr->color == RED) {
                                /*
                                 * Red node having red child. B- black, R-red
                                 *     B                R
                                 *    / \             /   \
                                 *   R   R  =>     B     B
                                 *  /               /   
                                 * R               R
                                 */
                                stack[ht - 2]->color = RED;
                                stack[ht - 1]->color = yPtr->color = BLACK;
                                ht = ht -2;
                        } else {
                                if (dir[ht - 1] == 0) {
                                        yPtr = stack[ht - 1];
                                } else {
                                        /*
                                         * XR - node X with red color
                                         * YR - node Y with red color
                                         * Red node having red child
                                         *(do single rotation left b/w X and Y)
                                         *         B             B
                                         *        /             /
                                         *      XR     =>      YR
                                         *        \           /
                                         *         YR        XR
                                         * one more additional processing will be
                                         * performed after this else part.  Since
                                         * we have red node (YR) with red child(XR)
                                         */
                                        xPtr = stack[ht - 1];
                                        yPtr = xPtr->link[1];
                                        xPtr->link[1] = yPtr->link[0];
                                        yPtr->link[0] = xPtr;
                                        stack[ht - 2]->link[0] = yPtr;
                                }
                                /*
                                 *  Red node(YR) with red child (XR) - single
                                 *  rotation b/w YR and XR for height balance. Still,
                                 *  red node (YR) is having red child.  So, change the
                                 *  color of Y to black and Black child B to Red R
                                 *          B           YR          YB
                                 *         /           /  \        /  \
                                 *        YR  =>   XR   B  =>  XR  R
                                 *       /
                                 *      XR
                                 */
                                xPtr = stack[ht - 2];
                                xPtr->color = RED;
                                yPtr->color = BLACK;
                                xPtr->link[0] = yPtr->link[1];
                                yPtr->link[1] = xPtr;
                                if (xPtr == root) {
                                        root = yPtr;
                                } else {
                                        stack[ht - 3]->link[dir[ht - 3]] = yPtr;
                                }
                                break;
                        }
                } else {
                        yPtr = stack[ht - 2]->link[0];
                        if ((yPtr != NULL) && (yPtr->color == RED)) {
                                /*
                                 * Red node with red child
                                 *        B             R
                                 *      /   \         /   \
                                 *     R     R =>  B     B
                                 *             \              \
                                 *              R              R
                                 *
                                 */
                                stack[ht - 2]->color = RED;
                                stack[ht - 1]->color = yPtr->color = BLACK;
                                ht = ht - 2;
                        } else {
                                if (dir[ht - 1] == 1) {
                                        yPtr = stack[ht - 1];
                                } else {
                                        /*
                                         * Red node(XR) with red child(YR) 
                                         *   B          B
                                         *    \          \
                                         *     XR  => YR
                                         *    /            \
                                         *   YR             XR
                                         * Single rotation b/w XR(node x with red color) & YR
                                         */
                                        xPtr = stack[ht - 1];
                                        yPtr = xPtr->link[0];
                                        xPtr->link[0] = yPtr->link[1];
                                        yPtr->link[1] = xPtr;
                                        stack[ht - 2]->link[1] = yPtr;
                                }
                                /*
                                 *   B              YR          YB
                                 *    \             /  \        /  \
                                 *     YR  =>   B   XR => R    XR
                                 *      \
                                 *       XR
                                 * Single rotation b/w YR and XR and change the color to
                                 * satisfy rebalance property.
                                 */
                                xPtr = stack[ht - 2];
                                yPtr->color = BLACK;
                                xPtr->color = RED;
                                xPtr->link[1] = yPtr->link[0];
                                yPtr->link[0] = xPtr;
                                if (xPtr == root) {
                                        root = yPtr;
                                } else {
                                        stack[ht - 3]->link[dir[ht - 3]] = yPtr;
                                }
                                break;
                        }
                }
        }
        root->color = BLACK;
}

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ConnectedVertNode tree  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//======================================================================================================================================
  
void fixup_ConnectedVertNode(struct VertexNode * Vptr,struct ConnectedVertNode *stack[98], int dir[98], int ht)
{

	struct ConnectedVertNode *xPtr, *yPtr;
        while ((ht >= 3) && (stack[ht - 1]->color == RED)) {
                if (dir[ht - 2] == 0) {
                        yPtr = stack[ht - 2]->link[1];
                        if (yPtr != NULL && yPtr->color == RED) {
                                /*
                                 * Red node having red child. B- black, R-red
                                 *     B                R
                                 *    / \             /   \
                                 *   R   R  =>     B     B
                                 *  /               /   
                                 * R               R
                                 */
                                stack[ht - 2]->color = RED;
                                stack[ht - 1]->color = yPtr->color = BLACK;
                                ht = ht -2;
                        } else {
                                if (dir[ht - 1] == 0) {
                                        yPtr = stack[ht - 1];
                                } else {
                                        /*
                                         * XR - node X with red color
                                         * YR - node Y with red color
                                         * Red node having red child
                                         *(do single rotation left b/w X and Y)
                                         *         B             B
                                         *        /             /
                                         *      XR     =>      YR
                                         *        \           /
                                         *         YR        XR
                                         * one more additional processing will be
                                         * performed after this else part.  Since
                                         * we have red node (YR) with red child(XR)
                                         */
                                        xPtr = stack[ht - 1];
                                        yPtr = xPtr->link[1];
                                        xPtr->link[1] = yPtr->link[0];
                                        yPtr->link[0] = xPtr;
                                        stack[ht - 2]->link[0] = yPtr;
                                }
                                /*
                                 *  Red node(YR) with red child (XR) - single
                                 *  rotation b/w YR and XR for height balance. Still,
                                 *  red node (YR) is having red child.  So, change the
                                 *  color of Y to black and Black child B to Red R
                                 *          B           YR          YB
                                 *         /           /  \        /  \
                                 *        YR  =>   XR   B  =>  XR  R
                                 *       /
                                 *      XR
                                 */
                                xPtr = stack[ht - 2];
                                xPtr->color = RED;
                                yPtr->color = BLACK;
                                xPtr->link[0] = yPtr->link[1];
                                yPtr->link[1] = xPtr;
                                if (xPtr == Vptr->root_cvert) {
                                        Vptr->root_cvert = yPtr;
                                } else {
                                        stack[ht - 3]->link[dir[ht - 3]] = yPtr;
                                }
                                break;
                        }
                } else {
                        yPtr = stack[ht - 2]->link[0];
                        if ((yPtr != NULL) && (yPtr->color == RED)) {
                                /*
                                 * Red node with red child
                                 *        B             R
                                 *      /   \         /   \
                                 *     R     R =>  B     B
                                 *             \              \
                                 *              R              R
                                 *
                                 */
                                stack[ht - 2]->color = RED;
                                stack[ht - 1]->color = yPtr->color = BLACK;
                                ht = ht - 2;
                        } else {
                                if (dir[ht - 1] == 1) {
                                        yPtr = stack[ht - 1];
                                } else {
                                        /*
                                         * Red node(XR) with red child(YR) 
                                         *   B          B
                                         *    \          \
                                         *     XR  => YR
                                         *    /            \
                                         *   YR             XR
                                         * Single rotation b/w XR(node x with red color) & YR
                                         */
                                        xPtr = stack[ht - 1];
                                        yPtr = xPtr->link[0];
                                        xPtr->link[0] = yPtr->link[1];
                                        yPtr->link[1] = xPtr;
                                        stack[ht - 2]->link[1] = yPtr;
                                }
                                /*
                                 *   B              YR          YB
                                 *    \             /  \        /  \
                                 *     YR  =>   B   XR => R    XR
                                 *      \
                                 *       XR
                                 * Single rotation b/w YR and XR and change the color to
                                 * satisfy rebalance property.
                                 */
                                xPtr = stack[ht - 2];
                                yPtr->color = BLACK;
                                xPtr->color = RED;
                                xPtr->link[1] = yPtr->link[0];
                                yPtr->link[0] = xPtr;
                                if (xPtr == Vptr->root_cvert) {
                                        Vptr->root_cvert = yPtr;
                                } else {
                                        stack[ht - 3]->link[dir[ht - 3]] = yPtr;
                                }
                                break;
                        }
                }
        }
        Vptr->root_cvert->color = BLACK;
  }
  
  
  struct ConnectedVertNode * createConnectedVertNode(int ver_no) 
  {
  	NE++; //whenever a newnode is created in conncetedVert tree in any Vertex node:NE++
        struct ConnectedVertNode *newnode;
        newnode = (struct ConnectedVertNode *)malloc(sizeof(struct ConnectedVertNode));
        newnode->conn_vert_no = ver_no;
        
        newnode->color = RED;
        newnode->link[0] = newnode->link[1] = NULL;
        return newnode;
  }
  
  
void insert_conn_vert(struct VertexNode * Vptr,int ver_no)
{ 
	//struct ConnectedVertNode * root_cvert1=Vptr->root_cvert;
	
  	//if(root_cvert!=NULL){printf("..Not Null..\n");}
        struct ConnectedVertNode *stack[98], *ptr, *newnode;
        int dir[98], ht = 0, index;
        int node_verno, node_tnc,node_vnc,tnc,vnc;
        double x,y,z,x1,y1,z1,x2,y2,z2;
        double lenE, ep_theta;
        
        x=Vptr->x; y=Vptr->y;z=Vptr->z; //(x,y,z) of the VertexNode itself. Needed to calculate the lenE of the edge
        
        tnc=ver_no/3; vnc=ver_no%3;
        x1=TV[tnc][vnc][0]; y1=TV[tnc][vnc][1]; z1=TV[tnc][vnc][2];
        
        ptr = Vptr->root_cvert;
        if (!Vptr->root_cvert) 
        {
        	//printf("Null.. ");
                Vptr->root_cvert = createConnectedVertNode(ver_no);
                Vptr->root_cvert->color = BLACK;
                
                return;
        }
        //printf("entered...\n");
        stack[ht] = Vptr->root_cvert;
        dir[ht++] = 0;
        /* find the place to insert the new node */
        while (ptr != NULL) 
        {
        	
        	node_verno=ptr->conn_vert_no; node_tnc=node_verno/3; node_vnc=node_verno%3;
        	//printf("enterd..ver_no=%d, %d, %d, node_verno=%d, %d, %d \n",ver_no,tnc,vnc,node_verno,node_tnc,node_vnc);
        	x2=TV[node_tnc][node_vnc][0]; y2=TV[node_tnc][node_vnc][1]; z2=TV[node_tnc][node_vnc][2];
        	
                if (x2==x1 && y2==y1 && z2==z1) 
                {
                
                 	//duplicate value not allowed 
  			NE1++;
  			lenE=sqrt((x2-x)*(x2-x)+(y2-y)*(y2-y)+(z2-z)*(z2-z));
  			ep_theta=get_ep_theta(tnc,node_tnc);
  			//printf("lenE=%lf, ep_theta=%lf \n",lenE,ep_theta);
  			imc2=imc2+0.5*lenE*ep_theta;
  			     
                        //return root_cvert;
                        return;
                }
                
                
                if(x2>x1) {index=1;}
                else if (x2<x1){index=0;}
                else
                {
                	if(y2>y1){index=1;}
                	else if(y2<y1){index=0;}
                	else
                	{
                		if(z2>z1){index=1;}
                		else if(z2<z1){index=0;}
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
        stack[ht - 1]->link[index] = newnode = createConnectedVertNode(ver_no); 
        fixup_ConnectedVertNode(Vptr,stack,dir,ht);
        //return root_cvert;
}

void free_ConnectedVertNode(struct ConnectedVertNode * start_root)
{
	if (start_root == NULL) {return;}
	free_ConnectedVertNode(start_root->link[0]);
	free_ConnectedVertNode(start_root->link[1]);
	
	free(start_root);
	

}

void free_VertexNode(struct VertexNode * start_root)
{
	if (start_root == NULL) {return;}
	free_VertexNode(start_root->link[0]);
	
	free_ConnectedVertNode(start_root->root_cvert);
	start_root->root_cvert=NULL;
	
	free_VertexNode(start_root->link[1]);
	
	free(start_root);


}
  


  
