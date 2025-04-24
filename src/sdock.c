/*                              SDOCK2.0
 *   Global Protein-Protein Docking Using Stepwise Force-Field Potentials   
 *                Written by ZHANG Changsheng 
 *  Copyright (c) 2024 Molecular Design Lab at Peking University, Beijing, China  
 
 *    contact: changshengzhang@pku.edu.cn, lhlai@pku.edu.cn
 *
 * sock.c : Map the two proteins to grids and do protein-protein global docking by FFT method.
 ***************************************************************************/


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
# include <time.h>
# include <sys/types.h>
# include <sys/times.h>
#include "mem.h"
#include "fftmem.h"
#include "grid.h"
#include "record.h"
#include "cluster.h"
#include "rotplan.h"

void current_cpu_time(int *h,int *m,int *s,int *ms)
{
          struct tms CPU_current;
          times(&CPU_current);
          *h = CPU_current.tms_utime / 360000;
          *m = CPU_current.tms_utime % 360000;
          *s = *m % 6000;
          *m = *m / 6000;
          *ms = *s % 100;
          *s = *s / 100;
}

fftwf_plan p_Rrepgrid,p_Rvdwgrid,p_Relegrid,p_Rsolgrid1,p_Rsolgrid2,p_Rsolgrid3,p_Rwatgrid,p_Rindgrid,p_RNhbgrid,p_ROhbgrid;
fftwf_plan p_Lrepgrid,p_Lvdwgrid,p_Lelegrid,p_Lsolgrid1,p_Lsolgrid2,p_Lsolgrid3,p_Lwatgrid,p_Lindgrid,p_LNhbgrid,p_LOhbgrid;
fftwf_plan p_ffttmp1,p_ffttmp;

void makeplan(int wmod)
{
	p_Rrepgrid=fftwf_plan_dft_r2c_3d(GRIDSIZE[0],GRIDSIZE[1],GRIDSIZE[2],Rrepgrid,Rrepgrid2,FFTW_ESTIMATE);
	p_Rvdwgrid=fftwf_plan_dft_3d(GRIDSIZE[0],GRIDSIZE[1],GRIDSIZE[2],Rvdwgrid,Rvdwgrid,FFTW_BACKWARD,FFTW_ESTIMATE);
	p_Relegrid=fftwf_plan_dft_3d(GRIDSIZE[0],GRIDSIZE[1],GRIDSIZE[2],Relegrid,Relegrid,FFTW_BACKWARD,FFTW_ESTIMATE);
	p_Rsolgrid1=fftwf_plan_dft_3d(GRIDSIZE[0],GRIDSIZE[1],GRIDSIZE[2],Rsolgrid1,Rsolgrid1,FFTW_BACKWARD,FFTW_ESTIMATE);
	p_Rsolgrid2=fftwf_plan_dft_3d(GRIDSIZE[0],GRIDSIZE[1],GRIDSIZE[2],Rsolgrid2,Rsolgrid2,FFTW_BACKWARD,FFTW_ESTIMATE);
	p_Rsolgrid3=fftwf_plan_dft_3d(GRIDSIZE[0],GRIDSIZE[1],GRIDSIZE[2],Rsolgrid3,Rsolgrid3,FFTW_BACKWARD,FFTW_ESTIMATE);
	if(wmod==WATMODE) p_Rwatgrid=fftwf_plan_dft_3d(GRIDSIZE[0],GRIDSIZE[1],GRIDSIZE[2],Rwatgrid,Rwatgrid,FFTW_BACKWARD,FFTW_ESTIMATE);
	p_Rindgrid=fftwf_plan_dft_3d(GRIDSIZE[0],GRIDSIZE[1],GRIDSIZE[2],Rindgrid,Rindgrid,FFTW_BACKWARD,FFTW_ESTIMATE);
	p_RNhbgrid=fftwf_plan_dft_3d(GRIDSIZE[0],GRIDSIZE[1],GRIDSIZE[2],RNhbgrid,RNhbgrid,FFTW_BACKWARD,FFTW_ESTIMATE);
	p_ROhbgrid=fftwf_plan_dft_3d(GRIDSIZE[0],GRIDSIZE[1],GRIDSIZE[2],ROhbgrid,ROhbgrid,FFTW_BACKWARD,FFTW_ESTIMATE);
	
	p_Lrepgrid=fftwf_plan_dft_r2c_3d(GRIDSIZE[0],GRIDSIZE[1],GRIDSIZE[2],Lrepgrid,Lrepgrid2,FFTW_MEASURE);
	p_Lvdwgrid=fftwf_plan_dft_3d(GRIDSIZE[0],GRIDSIZE[1],GRIDSIZE[2],Lvdwgrid,Lvdwgrid,FFTW_FORWARD,FFTW_MEASURE);
	p_Lelegrid=fftwf_plan_dft_3d(GRIDSIZE[0],GRIDSIZE[1],GRIDSIZE[2],Lelegrid,Lelegrid,FFTW_FORWARD,FFTW_MEASURE);
	p_Lsolgrid1=fftwf_plan_dft_3d(GRIDSIZE[0],GRIDSIZE[1],GRIDSIZE[2],Lsolgrid1,Lsolgrid1,FFTW_FORWARD,FFTW_MEASURE);
	p_Lsolgrid2=fftwf_plan_dft_3d(GRIDSIZE[0],GRIDSIZE[1],GRIDSIZE[2],Lsolgrid2,Lsolgrid2,FFTW_FORWARD,FFTW_MEASURE);
	p_Lsolgrid3=fftwf_plan_dft_3d(GRIDSIZE[0],GRIDSIZE[1],GRIDSIZE[2],Lsolgrid3,Lsolgrid3,FFTW_FORWARD,FFTW_MEASURE);
	if(wmod==WATMODE) p_Lwatgrid=fftwf_plan_dft_3d(GRIDSIZE[0],GRIDSIZE[1],GRIDSIZE[2],Lwatgrid,Lwatgrid,FFTW_FORWARD,FFTW_MEASURE);
	p_Lindgrid=fftwf_plan_dft_3d(GRIDSIZE[0],GRIDSIZE[1],GRIDSIZE[2],Lindgrid,Lindgrid,FFTW_FORWARD,FFTW_MEASURE);
	p_LNhbgrid=fftwf_plan_dft_3d(GRIDSIZE[0],GRIDSIZE[1],GRIDSIZE[2],LNhbgrid,LNhbgrid,FFTW_FORWARD,FFTW_MEASURE);
	p_LOhbgrid=fftwf_plan_dft_3d(GRIDSIZE[0],GRIDSIZE[1],GRIDSIZE[2],LOhbgrid,LOhbgrid,FFTW_FORWARD,FFTW_MEASURE);
	
	p_ffttmp1=fftwf_plan_dft_c2r_3d(GRIDSIZE[0],GRIDSIZE[1],GRIDSIZE[2],ffttmp,repscore,FFTW_MEASURE);
	p_ffttmp=fftwf_plan_dft_3d(GRIDSIZE[0],GRIDSIZE[1],GRIDSIZE[2],ffttmp,ffttmp,FFTW_BACKWARD,FFTW_MEASURE);
}
void detroyplan(int wmod)
{
  	fftwf_destroy_plan( p_Rrepgrid) ;
  	fftwf_destroy_plan( p_Rvdwgrid ) ;
   	fftwf_destroy_plan( p_Relegrid ) ;
 	fftwf_destroy_plan( p_Rsolgrid1 ) ;
  	fftwf_destroy_plan( p_Rsolgrid2) ;
  	fftwf_destroy_plan( p_Rsolgrid3) ;
  	if(wmod==WATMODE) fftwf_destroy_plan( p_Rwatgrid) ;
  	fftwf_destroy_plan( p_Rindgrid) ;
  	fftwf_destroy_plan( p_RNhbgrid) ;
  	fftwf_destroy_plan( p_ROhbgrid) ;

   	fftwf_destroy_plan( p_Lrepgrid) ;
 	fftwf_destroy_plan( p_Lvdwgrid ) ;
    	fftwf_destroy_plan( p_Lelegrid ) ;
 	fftwf_destroy_plan( p_Lsolgrid1 ) ;
  	fftwf_destroy_plan( p_Lsolgrid2) ;
  	fftwf_destroy_plan( p_Lsolgrid3) ;
  	if(wmod==WATMODE) fftwf_destroy_plan( p_Lwatgrid) ;
  	fftwf_destroy_plan( p_Lindgrid) ;
  	fftwf_destroy_plan( p_LNhbgrid) ;
  	fftwf_destroy_plan( p_LOhbgrid) ;
	
  	fftwf_destroy_plan( p_ffttmp1 ) ;
  	fftwf_destroy_plan( p_ffttmp) ;	
}
void receptorFFT(int wmod)
{
	fftwf_execute(p_Rrepgrid) ;
	fftwf_execute(p_Rvdwgrid) ;
	fftwf_execute(p_Relegrid) ;
 	fftwf_execute(p_Rsolgrid1 ) ;
  	fftwf_execute(p_Rsolgrid2) ;
	fftwf_execute(p_Rsolgrid3) ;
	if(wmod==WATMODE) fftwf_execute(p_Rwatgrid) ;
	fftwf_execute(p_Rindgrid) ;
	fftwf_execute(p_RNhbgrid) ;
	fftwf_execute(p_ROhbgrid) ;
}
void ligandFFT(int wmod)
{
	fftwf_execute(p_Lrepgrid) ;
	fftwf_execute(p_Lvdwgrid) ;
	fftwf_execute(p_Lelegrid) ;
 	fftwf_execute(p_Lsolgrid1 ) ;
  	fftwf_execute(p_Lsolgrid2) ;
	fftwf_execute(p_Lsolgrid3) ;
	if(wmod==WATMODE) fftwf_execute(p_Lwatgrid) ;
	fftwf_execute(p_Lindgrid) ;
	fftwf_execute(p_LNhbgrid) ;
	fftwf_execute(p_LOhbgrid) ;
}
void calFFTmult(fftwf_complex *fg1,fftwf_complex *fg2,fftwf_complex *fg3)
{
    int i;
          
	for(i=0;i<GRID3SIZE;i++){
		fg3[i][0]=fg1[i][0]*fg2[i][0]-fg1[i][1]*fg2[i][1];
		fg3[i][1]=fg1[i][1]*fg2[i][0]+fg1[i][0]*fg2[i][1] ;
	}
}
void calFFTmult2(fftwf_complex *fg1,fftwf_complex *fg2,fftwf_complex *fg3)
{
    int i;
          
	for(i=0;i<GRID3SIZE;i++){
		fg3[i][0]=fg1[i][0]*fg2[i][0]+fg1[i][1]*fg2[i][1];
		fg3[i][1]=fg1[i][1]*fg2[i][0]-fg1[i][0]*fg2[i][1] ;
	}
}
void inittot()
{
	int i,j,k,g1;
	int ii,jj,kk,g2;
	
	for(i=0;i<GRIDSIZE[0];i++){
		if(i!=0) ii=GRIDSIZE[0]-i;
		else ii=0;
		for(j=0;j<GRIDSIZE[1];j++){
			if(j!=0) jj=GRIDSIZE[1]-j;
			else jj=0;
			for(k=0;k<GRIDSIZE[2];k++){
				if(k!=0) kk=GRIDSIZE[2]-k;
				else kk=0;
				g1=getgindex(i,j,k);
				g2=getgindex(ii,jj,kk);
				totalscore[g2]=repscore[g1]*REPWEIGHT/(GRID3SIZE);
			}
		}
	}
}
void inverseFFT(int wmod)
{
	int i;
      float tden;

      tden=3.0/4.0*GRIDSPACESIZE*GRIDSPACESIZE*GRIDSPACESIZE/(SURFWATRAD*SURFWATRAD*SURFWATRAD*3.14159265358979323846);

	calFFTmult2(Rrepgrid2,Lrepgrid2,ffttmp);
	fftwf_execute(p_ffttmp1) ;
	inittot();
	for(i=0;i<GRID3SIZE;i++)
		repscore[i]=totalscore[i]/REPWEIGHT;
	calFFTmult(Rvdwgrid,Lvdwgrid,ffttmp);
	fftwf_execute(p_ffttmp);
	for(i=0;i<GRID3SIZE;i++){
		vdwscore[i]=ffttmp[i][1]/(GRID3SIZE)/2.0;
		totalscore[i]+=vdwscore[i];
	}
	calFFTmult(Relegrid,Lelegrid,ffttmp);
	fftwf_execute(p_ffttmp);
	for(i=0;i<GRID3SIZE;i++){
		elescore[i]=ffttmp[i][1]/(GRID3SIZE)/2.0;
		totalscore[i]+=elescore[i]*ELEWEIGHT;
	}
	calFFTmult(Rsolgrid1,Lsolgrid1,ffttmp);
	fftwf_execute(p_ffttmp);
	for(i=0;i<GRID3SIZE;i++){
		solscore[i]=ffttmp[i][1]/(GRID3SIZE);
	}
	calFFTmult(Rsolgrid2,Lsolgrid2,ffttmp);
	fftwf_execute(p_ffttmp);
	for(i=0;i<GRID3SIZE;i++){
		solscore[i]+=ffttmp[i][1]/(GRID3SIZE);
	}
	calFFTmult(Rsolgrid3,Lsolgrid3,ffttmp);
	fftwf_execute(p_ffttmp);
	for(i=0;i<GRID3SIZE;i++){
		solscore[i]+=ffttmp[i][1]/(GRID3SIZE);
		solscore[i]=solscore[i]/2.0;
		totalscore[i]+=solscore[i]*SOLWEIGHT;
	}
	if(wmod==WATMODE) {
		calFFTmult(Rwatgrid,Lwatgrid,ffttmp);
		fftwf_execute(p_ffttmp);
		for(i=0;i<GRID3SIZE;i++){
			watscore[i]=-(ffttmp[i][1]/(GRID3SIZE)*tden);
			totalscore[i]+=watscore[i]*WATWEIGHT;
		}
	}
	calFFTmult(Rindgrid,Lindgrid,ffttmp);
	fftwf_execute(p_ffttmp);
	for(i=0;i<GRID3SIZE;i++){
		indscore[i]=ffttmp[i][1]/(GRID3SIZE);
		totalscore[i]+=indscore[i]*INDWEIGHT;
	}
	calFFTmult(RNhbgrid,LNhbgrid,ffttmp);
	fftwf_execute(p_ffttmp);
	for(i=0;i<GRID3SIZE;i++){
		Nhbscore[i]=ffttmp[i][1]/(GRID3SIZE);
		totalscore[i]+=Nhbscore[i]*HBNWEIGHT;
	}
	calFFTmult(ROhbgrid,LOhbgrid,ffttmp);
	fftwf_execute(p_ffttmp);
	for(i=0;i<GRID3SIZE;i++){
		Ohbscore[i]=ffttmp[i][1]/(GRID3SIZE);
		totalscore[i]+=Ohbscore[i]*HBOWEIGHT;
	}
}
#define KEEPSCORE  -8.0
void record_best(int roti)
{
	int x,y,z,xyz;
	int fx,fy,fz;

	for( x = 0 ; x < GRIDSIZE[0] ; x ++ ) {
		fx = x ;
		if( fx > ( GRIDSIZE[0] / 2 ) ) fx -= GRIDSIZE[0] ;

		for( y = 0 ; y < GRIDSIZE[1] ; y ++ ) {
			fy = y ;
			if( fy > ( GRIDSIZE[1] / 2 ) ) fy -= GRIDSIZE[1] ;

			for( z = 0 ; z < GRIDSIZE[2] ; z ++ ) {
				fz = z ;
				if( fz > ( GRIDSIZE[2] / 2 ) ) fz -= GRIDSIZE[2] ;

				xyz = z + GRIDSIZE[2]  * ( y + GRIDSIZE[1] * x ) ;
				if(totalscore[xyz]<KEEPSCORE){
					save_dock(roti, fx*GRIDSPACESIZE,fy*GRIDSPACESIZE,fz*GRIDSPACESIZE,totalscore[xyz],repscore[xyz],vdwscore[xyz],elescore[xyz],solscore[xyz],watscore[xyz],indscore[xyz],Nhbscore[xyz],Ohbscore[xyz]);
				}
			}
		}
	}
}

void Dock(STRUCTURE *ligand,SURFWAT *ligwat,float prog_per,int wmod)
{
	int ri;
	float m[3][3];
	int t=(ROTN*prog_per*0.01+1);

	for(ri=0;ri<ROTN;ri++){
		if(ri%t==0||ri==ROTN-1)
			printf("%6.2f\n",(ri+1)*100.0/ROTN);
		quat2matrix(ROTPLAN[ri], m);
		rotligand(ligand,m);
		if(wmod==WATMODE) rotligwat(ligwat,m);
		gen_ligand_grid(ligand,ligwat,wmod);
		ligandFFT(wmod);
		inverseFFT(wmod);
		record_best(ri);
	}
}

void create_grid(int wmod)
{
	fftnew(mattergrid,GRID3SIZE);
	fftnew(Rrepgrid,GRID3SIZE);
	fftnew(Rrepgrid2,GRID3SIZE);
	fftnew(Rvdwgrid,GRID3SIZE);
	fftnew(Relegrid,GRID3SIZE);
	fftnew(Rsolgrid1,GRID3SIZE);
	fftnew(Rsolgrid2,GRID3SIZE);
	fftnew(Rsolgrid3,GRID3SIZE);
	if(wmod==WATMODE) fftnew(Rwatgrid,GRID3SIZE);
	fftnew(Rindgrid,GRID3SIZE);
	fftnew(RNhbgrid,GRID3SIZE);
	fftnew(ROhbgrid,GRID3SIZE);
	fftnew(Lrepgrid,GRID3SIZE);
	fftnew(Lrepgrid2,GRID3SIZE);
	fftnew(Lvdwgrid,GRID3SIZE);
	fftnew(Lelegrid,GRID3SIZE);
	fftnew(Lsolgrid1,GRID3SIZE);
	fftnew(Lsolgrid2,GRID3SIZE);
	fftnew(Lsolgrid3,GRID3SIZE);
	if(wmod==WATMODE) fftnew(Lwatgrid,GRID3SIZE);
	fftnew(Lindgrid,GRID3SIZE);
	fftnew(LNhbgrid,GRID3SIZE);
	fftnew(LOhbgrid,GRID3SIZE);
	fftnew(repscore,GRID3SIZE);
	fftnew(vdwscore,GRID3SIZE);
	fftnew(elescore,GRID3SIZE);
	fftnew(solscore,GRID3SIZE);
	fftnew(watscore,GRID3SIZE);
	fftnew(indscore,GRID3SIZE);
	fftnew(Nhbscore,GRID3SIZE);
	fftnew(Ohbscore,GRID3SIZE);
	fftnew(ffttmp,GRID3SIZE);
	snew(totalscore,GRID3SIZE);
}
void destroy_grid(int wmod)
{
	fftfree(mattergrid);
	fftfree(Rrepgrid);
	fftfree(Rrepgrid2);
	fftfree(Rvdwgrid);
	fftfree(Relegrid);
	fftfree(Rsolgrid1);
	fftfree(Rsolgrid2);
	fftfree(Rsolgrid3);
	if(wmod==WATMODE) fftfree(Rwatgrid);
	fftfree(Rindgrid);
	fftfree(RNhbgrid);
	fftfree(ROhbgrid);
	fftfree(Lrepgrid);
	fftfree(Lrepgrid2);
	fftfree(Lvdwgrid);
	fftfree(Lelegrid);
	fftfree(Lsolgrid1);
	fftfree(Lsolgrid2);
	fftfree(Lsolgrid3);
	if(wmod==WATMODE) fftfree(Lwatgrid);
	fftfree(Lindgrid);
	fftfree(LNhbgrid);
	fftfree(LOhbgrid);
	fftfree(repscore);
	fftfree(vdwscore);
	fftfree(elescore);
	fftfree(solscore);
	fftfree(watscore);
	fftfree(indscore);
	fftfree(Nhbscore);
	fftfree(Ohbscore);
	fftfree(ffttmp);
	sfree(totalscore);
}

void print_sdock_usage()
{
	printf(" Usage: sdock preprocessed_protein_structure1 preprocessed_protein_structure2 \n");
	printf("    [watermap_for_protein1] [watermap_for_protein2] \n");
	printf("    [-o sdock_result_file_name (sdock_record)]  \n");
	printf("    [-p progress with percentage (0.1)]  \n");
	printf("    [-r rotation_sampling_file (so3layer.qua)] \n");
	printf("    [-c weight_of_collision_term (0.17)]  \n");
	printf("    [-e weight_of_electrostatic_term (1.4)]\n");
	printf("    [-s weight_of_desolvation_term (0.67/0.60)]  \n");
	printf("    [-i weight_of_induction_term (4.0)]  \n");
	printf("    [-H weight_of_backboneO_HB_term] (7.0)\n");
	printf("    [-h weight_of_backboneN_HB_term (0.5/1.0)]  \n");
	printf("    [-w weight_of_explicit_water_term (0.0/1.0)]  \n");
	printf("    [-n cluster_number (1000)]  \n");
	printf("    [-d cluster_radius (3.0)]\n\n");
	printf(" Examples:\n");
	printf(" sdock 1AY7R.pdb 1AY7L.pdb                       \n      #no watmap, explicit water term is not included\n");
	printf(" sdock 1AY7R.pdb 1AY7L.pdb 1AY7Rwat.pdb 1AY7Lwat.pdb\n      #watmaps are inputed, explicit water mode\n");
	printf(" sdock 1AY7R.pdb 1AY7L.pdb -o 1AY7_dock_result -p 0\n      #result file name and no progress is displayed\n");
	printf(" sdock 1AY7R.pdb 1AY7L.pdb -o 1AY7_dock_result -H 8.0\n      #a larger weight for backbone oxygen hydrogen bond term\n");
	printf(" sdock 1AY7R.pdb 1AY7L.pdb 1AY7Rwat.pdb 1AY7Lwat.pdb -s 0.63\n      #a larger desolvation weight in explicit water mode\n");
	printf(" sdock 1AY7R.pdb 1AY7L.pdb -n 100 -d 2.0    \n      #only output the top 100 solutions with a smaller cluster radius\n");
}

int main(int argc, char *argv[])
{
	STRUCTURE proA,proB,*receptor,*ligand;
	SURFWAT watA,watB,*recwat=NULL,*ligwat=NULL;
	float prog_per;
	int arg_rotfile;
	float (*receptorca)[3],(*ligandca)[3];
	int receptorcaN,ligandcaN;
	int ih,im,is,ims;
	int t=0;
	int arg_i,arg_result;
	int dockmod=NOWATMODE;

	printf("                             SDOCK2.0\n");
	printf("  Global Protein-Protein Docking Using Stepwise Force-Field Potentials   \n");
	printf(" Copyright (c) 2024 Molecular design lab in Peking university, Beijing, China\n\n");
	printf(" sdock: protein-protein docking and save results.\n");
	
	if(argc==1){
		print_sdock_usage();
		return EXIT_SUCCESS;
	}

	if(argv[1][0]=='-'){
		printf("\nERROR: The first two options should be two structure files, and not %s!\n\n",argv[1]);
		print_sdock_usage();
		return EXIT_SUCCESS;
	}
	if(argc<=2){
		printf("\nERROR: You only input one structure file and you should input two!\n\n");
		print_sdock_usage();
		return EXIT_SUCCESS;
	}

	if(argv[2][0]=='-'){
		printf("\nERROR: The first two options should be two structure files, and not %s!\n\n",argv[2]);
		print_sdock_usage();
		return EXIT_SUCCESS;
	}
	if(argv[3][0]!='-'){
		dockmod=WATMODE;
		if(argv[4][0]=='-'){
			printf("\nERROR: You only input one watermap file and you should input two!\n\n");
			print_sdock_usage();
			return EXIT_SUCCESS;
		}
	}

	for(arg_i=3+dockmod*2;arg_i<argc;arg_i+=2){
		if(argv[arg_i][0]!='-'||(argv[arg_i][1]!='o'&&argv[arg_i][1]!='p'&&argv[arg_i][1]!='r'&&argv[arg_i][1]!='c'&&argv[arg_i][1]!='e'&&argv[arg_i][1]!='s'&&argv[arg_i][1]!='i'&&argv[arg_i][1]!='H'&&argv[arg_i][1]!='h'&&argv[arg_i][1]!='w'&&argv[arg_i][1]!='n'&&argv[arg_i][1]!='d')){
			printf("\nERROR at '%s': should be one of these options: -o -p -r -c -e -s -i -H -h -w -n -d!\n\n",argv[arg_i]);
			print_sdock_usage();
			return EXIT_SUCCESS;
		}
		if(arg_i+1>=argc||argv[arg_i+1][0]=='-'){
			printf("\nplease input the option value after %s!\n\n",argv[arg_i]);
			print_sdock_usage();
			return EXIT_SUCCESS;
		}
	}

	arg_result=0;
	prog_per=0.1;
	arg_rotfile=0;
      REPWEIGHT=0.17;
      ELEWEIGHT=1.4;
      if(dockmod==NOWATMODE)SOLWEIGHT=0.67; else SOLWEIGHT=0.60;
      if(dockmod==NOWATMODE) WATWEIGHT=0.0; else WATWEIGHT=1.0;
      INDWEIGHT=4.0;
      HBOWEIGHT=7.0;
      if(dockmod==NOWATMODE) HBNWEIGHT=0.5; else HBNWEIGHT=1.0;
	max_clusterN=1000;
	LRMSDCUTOFF=3.0;

	for(arg_i=3+dockmod*2;arg_i<argc;arg_i+=2){
		if(argv[arg_i][1]=='o')
			arg_result=arg_i+1;
		else if(argv[arg_i][1]=='p')
			prog_per=atof(argv[arg_i+1]);
		else if(argv[arg_i][1]=='r'){
			arg_rotfile=arg_i+1;
		}
		else if(argv[arg_i][1]=='c'){
			REPWEIGHT=atof(argv[arg_i+1]);
			if(REPWEIGHT>0.5){
				printf("ERROR: The weight of collision term (%5.3f) is too large!\n",REPWEIGHT);
				return EXIT_SUCCESS;
			}
			if(REPWEIGHT<0.00001)
				printf("WARNING: The weight of collision term cannot set to zero!\n");
		}
		else if(argv[arg_i][1]=='e'){
			ELEWEIGHT=atof(argv[arg_i+1]);
			if(ELEWEIGHT>5.0){
				printf("ERROR: The weight of electrostatic term (%5.3f) is too large!\n",ELEWEIGHT);
				return EXIT_SUCCESS;
			}
			if(ELEWEIGHT<0.00001)
				printf("WARNING: The weight of electrostatic term is set to zero!\n");
		}
		else if(argv[arg_i][1]=='s'){
			SOLWEIGHT=atof(argv[arg_i+1]);
			if(SOLWEIGHT>2.0){
				printf("ERROR: The weight of desolvation term (%5.3f) is too large!\n",SOLWEIGHT);
				return EXIT_SUCCESS;
			}
			if(SOLWEIGHT<0.00001)
				printf("WARNING: The weight of desolvation term is set to zero!\n");
		}
		else if(argv[arg_i][1]=='i'){
			INDWEIGHT=atof(argv[arg_i+1]);
			if(INDWEIGHT>10.0){
				printf("ERROR: The weight of induction term (%5.3f) is too large!\n",ELEWEIGHT);
				return EXIT_SUCCESS;
			}
			if(INDWEIGHT<0.00001)
				printf("WARNING: The weight of induction term is set to zero!\n");
		}
		else if(argv[arg_i][1]=='w'){
			WATWEIGHT=atof(argv[arg_i+1]);
			if(dockmod==WATMODE){
				if(WATWEIGHT>4.0){
					printf("ERROR: The weight of explicit water term (%5.3f) is too large!\n",WATWEIGHT);
					return EXIT_SUCCESS;
				}
				if(WATWEIGHT<0.00001)
					printf("WARNING: The weight of explicit water term is set to zero!\n");
			}
		}
		else if(argv[arg_i][1]=='H'){
			HBOWEIGHT=atof(argv[arg_i+1]);
			if(HBOWEIGHT>20.0){
				printf("ERROR: The weight of backboneO HB term (%5.3f) is too large!\n",HBOWEIGHT);
				return EXIT_SUCCESS;
			}
			if(HBOWEIGHT<0.00001)
				printf("WARNING: The weight of backboneO HB term is set to zero!\n");
		}
		else if(argv[arg_i][1]=='h'){
			HBNWEIGHT=atof(argv[arg_i+1]);
			if(HBNWEIGHT>20.0){
				printf("ERROR: The weight of backboneN HB term (%5.3f) is too large!\n",HBNWEIGHT);
				return EXIT_SUCCESS;
			}
			if(HBNWEIGHT<0.00001)
				printf("WARNING: The weight of backboneN HB term is set to zero!\n");
		}
		else if(argv[arg_i][1]=='n'){
			max_clusterN=atoi(argv[arg_i+1]);
			if(max_clusterN>9999){
				printf("WARNING: The cluster number (%d) is too large! and will be changed to 9999\n",max_clusterN);
				max_clusterN=9999;
			}
			if(max_clusterN==0){
				printf("ERROR: The cluster number can not be zero!\n");
				return EXIT_SUCCESS;
			}
		}
		else if(argv[arg_i][1]=='d'){
			LRMSDCUTOFF=atof(argv[arg_i+1]);
			if(LRMSDCUTOFF>10.0){
				printf("ERROR: The cluster radius (%5.3f) is too large!\n",LRMSDCUTOFF);
				return EXIT_SUCCESS;
			}
			if(LRMSDCUTOFF<0.00001)
				printf("WARNING: The cluster radius is set to zero!\n");
		}
	}

	if(LRMSDCUTOFF>=3.0)
		RESULTKEEPN=20*max_clusterN;
	else if(LRMSDCUTOFF>=2.5)
		RESULTKEEPN=12*max_clusterN;
	else if(LRMSDCUTOFF>=2.0)
		RESULTKEEPN=8*max_clusterN;
	else if(LRMSDCUTOFF>=1.5)
		RESULTKEEPN=4*max_clusterN;
	else
		RESULTKEEPN=2*max_clusterN;
	if(arg_rotfile==0)
		ROTN=readrotplan("so3layer.qua");
	else
		ROTN=readrotplan(argv[arg_rotfile]);


	init_record();
	read_structure(&proA,argv[1]);
	read_structure(&proB,argv[2]);
	if(dockmod==WATMODE){
		read_surfacewater(&watA, argv[3]);
		read_surfacewater(&watB, argv[4]);
	}

	t=decide_gridsize(&proA,&proB);
	if(t==0) {
		receptor = &proA;
		ligand = &proB;
		fixpro(receptor);
		if(dockmod==WATMODE){	
			recwat = &watA;
			ligwat = &watB;
			fixwat(recwat);	
		}
	}
	else{
		receptor = &proB;
		ligand = &proA;
		fixpro(receptor);	
		if(dockmod==WATMODE){	
			recwat = &watB;
			ligwat = &watA;
			fixwat(recwat);	
		}
	}

	create_grid(dockmod);

	current_cpu_time(&ih,&im,&is,&ims);
	makeplan(dockmod);
	gen_receptor_grid(receptor,recwat,dockmod);
	receptorFFT(dockmod);
	Dock(ligand,ligwat,prog_per,dockmod);
	detroyplan(dockmod);
	destroy_grid(dockmod);
	if(dockmod==WATMODE){
		sfree(recwat->w);
		sfree(ligwat->w);
	}

	snew(ligandca,ligand->CaN);	
	ligandcaN=get_ligandca(ligand, ligandca);
	snew(receptorca,receptor->CaN);	
	receptorcaN=get_ligandca(receptor, receptorca);

	arrange_record();
	cluster_result(receptorcaN,receptorca,ligandcaN,ligandca);

	if(arg_result==0&&arg_rotfile==0)
		print_cluster(dockmod,argv[1],argv[2],t,"sdock_record","so3layer.qua",REPWEIGHT,ELEWEIGHT,SOLWEIGHT,WATWEIGHT,INDWEIGHT,HBNWEIGHT,HBOWEIGHT);
	else if(arg_result==0&&arg_rotfile!=0)
		print_cluster(dockmod,argv[1],argv[2],t,"sdock_record",argv[arg_rotfile],REPWEIGHT,ELEWEIGHT,SOLWEIGHT,WATWEIGHT,INDWEIGHT,HBNWEIGHT,HBOWEIGHT);
	else if(arg_result!=0&&arg_rotfile==0)
		print_cluster(dockmod,argv[1],argv[2],t,argv[arg_result],"so3layer.qua",REPWEIGHT,ELEWEIGHT,SOLWEIGHT,WATWEIGHT,INDWEIGHT,HBNWEIGHT,HBOWEIGHT);
	else
		print_cluster(dockmod,argv[1],argv[2],t,argv[arg_result],argv[arg_rotfile],REPWEIGHT,ELEWEIGHT,SOLWEIGHT,WATWEIGHT,INDWEIGHT,HBNWEIGHT,HBOWEIGHT);


	sfree(receptorca);	
	sfree(ligandca);
	free_record();
	sfree(receptor->a);
	sfree(ligand->a);
	free_rotplan();
	current_cpu_time(&ih,&im,&is,&ims);
	printf("==> %s + %s DOCKING RUNTIME: %2d:%2d:%2d:%2d \n",argv[1],argv[2],ih,im,is,ims);

	return EXIT_SUCCESS;
}

