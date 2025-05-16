/*                              SDOCK2.0
 *   Global Protein-Protein Docking Using Stepwise Force-Field Potentials   
 *                Written by ZHANG Changsheng 
 *  Copyright (c) 2024 Molecular Design Lab at Peking University, Beijing, China  
 
 *    contact: changshengzhang@pku.edu.cn, lhlai@pku.edu.cn
 *
 * rotplan.h : read the rotation sampling file and translate the quaternion to rotational matrix.
 ***************************************************************************/
#ifndef _rotplan_h
#define _rotplan_h

int getrotplanNum(char *fn);
int readrotplan(char *fn,float (*ROTPLAN)[4]);
void quat2matrix(float quat[4], float m[3][3]);

#endif

