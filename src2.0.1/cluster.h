/*                              SDOCK2.0
 *   Global Protein-Protein Docking Using Stepwise Force-Field Potentials   
 *                Written by ZHANG Changsheng 
 *  Copyright (c) 2024 Molecular Design Lab at Peking University, Beijing, China  
 
 *    contact: changshengzhang@pku.edu.cn, lhlai@pku.edu.cn
 *
 * cluster.h : cluster the results by interface Calpha RMSD and output the docking results.
 ***************************************************************************/
#ifndef _cluster_h
#define _cluster_h

#include "record.h"

#define MAXMEMBER  300
#define MAXCLUSTER  1999
typedef struct{
	int N;
	float rmsdcutoff;
	int id[MAXMEMBER];
	float Lrmsd[MAXMEMBER];
}CLUSTER;

int cluster_result(int RCaN,float RCa[][3],int CaN,float Ca[][3], int max_clusterN, int RESULTKEEPN,CLUSTER cluster[MAXCLUSTER],DOCKRECORD *dock_result,int *result_index,float LRMSDCUTOFF,float (*ROTPLAN)[4]);
void print_cluster(int mode,char *pA, char *pB, int t, char *fn,char *rfn, float cw, float ew, float sw, float ww, float iw, float nw, float ow, int clusterN, CLUSTER cluster[MAXCLUSTER],DOCKRECORD *dock_result,int *result_index,float LRMSDCUTOFF);


#endif

