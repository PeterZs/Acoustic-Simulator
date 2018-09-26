#ifndef _CONSTLIST_H_
#define _CONSTLIST_H_

#define EDGE_NODE_NUM 		 2	
#define FACE_NODE_NUM     	 3
#define VOXEL_NODE_NUM    	 4

#define TOLERANCE 		  	 1e-5
#define PI 				  	 3.1415926535
#define INF 			  	 1e10

#define ERR_FILE_OPEN 	      1
#define ERR_FILE_FORMAT       2
#define ERR_WRONG_ARG_NUMBER  3

#define AREA_NORMALIZE 		  0
#define UNIT_NORMALIZE 		  1
#define REAL 				  3
#define IMG 				  4

#define READ_NODES 			  1
#define READ_VOXELS           2

enum ElementType{INNER, HARD, INLET, OUTLET};
enum Operation{MUFFLER_SIMULATION, EXPERIMENT};

#endif 