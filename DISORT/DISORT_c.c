/*
**	$Id: simple_c2f.c,v 1.2 1995/08/02 17:01:36 kirk Exp $
**
** NAME:
** 	disort_c.c
**
** PURPOSE:
**	This C function is used to demonstrate how to pass all IDL
**	simple varable types to a FORTRAN routine via a C wrapper function.
**	The passed variables are squared in the Fortran routine and
**	returned.
**
** CATEGORY:
**	Dynamic Link
**
*/

#include <stdio.h>

;/* Fortran routine */
extern "C"
{
  void disort_(int *, float *, float *, int *, float *, float *,
	       float *, float *, short *, int *, float *, int *,
	       short *, int *, float *, int *, float *, int *, float *,
	       float *, float *, float *, short *, float *, float *,
	       float *, float *, short *, short *, float *, short *,
	       short *, int *, int *, int *, int *,
	       int *, float *, float *, float *, float *, float *, float *,
	       float *, float *);

}

extern "C"
{
  void disort_c (int argc, void *argv[]);
}

void disort_c (int argc, void *argv[])
{

  /*
  ** Declare variables
  */
  short *USRTAU, *USRANG, *LAMBER, *PLANK, *ONLYFL, *PRNT, *HEADER;         /* Pointers to booleans       */ 
  int *NLYR, *NMOM, *NTAU, *NSTR, *NUMU, *NPHI, *IBCND;  
  int *MAXCLY, *MAXULV, *MAXUMU, *MAXPHI, *MAXMOM;        /* Pointers to integer       */  
  float *DTAUC, *SSALB, *PMOM, *TEMPER, *WVNMLO, *WVNMHI, *UTAU, *UMU;
  float *PHI, *FBEAM, *UMU0, *PHI0, *FISOT, *ALBEDO, *BTEMP, *TTEMP;
  float *TEMIS, *ACCUR, *RFLDIR, *RFLDN, * FLUP, *DFDT, *UAVG, *UU;
  float *ALBMED, *TRNMED;   
  ;                               /* Pointers to float       */

  /*
  ** Insure that the correct number of arguments were passed in (argc = 30)
  */

  if(argc != 45)
    {
      /*
      ** Print an error message and return
      */
      fprintf(stderr,"disort: Incorrect number of arguments\n");
    }
  /*
  ** Cast the pointer in argv to the pointer variables
  */
  NLYR  = (int *)   argv[0];
  DTAUC  = (float *)   argv[1];
  SSALB = (float *)   argv[2];
  NMOM = (int *)   argv[3];
  PMOM = (float *)   argv[4];
  TEMPER = (float *)   argv[5];
  WVNMLO  = (float *)   argv[6];
  WVNMHI  = (float *)   argv[7];
  USRTAU  = (short *)    argv[8];
  NTAU  = (int *)    argv[9];
  UTAU = (float *)    argv[10];
  NSTR  = (int *)    argv[11];
  USRANG  = (short *)    argv[12];
  NUMU = (int *)   argv[13];
  UMU = (float *)   argv[14];
  NPHI = (int *)   argv[15];
  PHI = (float *)   argv[16];
  IBCND = (int *)   argv[17];
  FBEAM = (float *)   argv[18];
  UMU0 = (float *)   argv[19];
  PHI0 = (float *)   argv[20];
  FISOT = (float *)   argv[21];
  LAMBER = (short *)   argv[22];
  ALBEDO = (float *)   argv[23];
  BTEMP = (float *)   argv[24];
  TTEMP = (float *)   argv[25];
  TEMIS = (float *)   argv[26];
  PLANK = (short *)   argv[27];
  ONLYFL = (short *)   argv[28];
  ACCUR = (float *)   argv[29];
  PRNT  = (short *)   argv[30];
  HEADER= (short *)   argv[31];
  MAXCLY= (int *)   argv[32];
  MAXULV= (int *)   argv[33];
  MAXUMU= (int *)   argv[34];
  MAXPHI= (int *)   argv[35];
  MAXMOM= (int *)   argv[36];
  RFLDIR= (float *)   argv[37];
  RFLDN= (float *)   argv[38];
  FLUP= (float *)   argv[39];
  DFDT= (float *)   argv[40];
  UAVG= (float *)   argv[41];
  UU= (float *)   argv[42];
  ALBMED= (float *)   argv[43];
  TRNMED= (float *)   argv[44];

  /*
  ** Now we need to call the fortran subroutine. The FORTRAN subroutine 
  ** uses varable length strings for the string parameter( CHARACTER*(*) ).
  ** Because of this we must pass in the length of each string that is 
  ** passed to the FORTRAN subprocedure. The string lengths are added 
  ** to the end of the parameter list. 
  */

  disort_(NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER,
	  WVNMLO, WVNMHI, USRTAU, NTAU, UTAU, NSTR,
	  USRANG, NUMU, UMU, NPHI, PHI, IBCND, FBEAM,
	  UMU0, PHI0, FISOT, LAMBER, ALBEDO, BTEMP,
	  TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, PRNT,
	  HEADER, MAXCLY, MAXULV, MAXUMU, MAXPHI,
	  MAXMOM, RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU,
	  ALBMED, TRNMED );

  /*
  ** That should be it
  */
} 

