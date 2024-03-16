#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<fftw3.h>
#include<omp.h>

#include"nbody.h"

/*------------------------------GLOBAL VARIABLES-------------------------------*/
//-----------------------------------------------------------------------------//
//                              from N-body code 
//-----------------------------------------------------------------------------//

float  vhh, // Hubble parameter in units of 100 km/s/Mpc
  vomegam, // Omega_matter, total matter density (baryons+CDM) parameter
  vomegalam, // Cosmological Constant 
  vomegab, //Omega_baryon
  sigma_8_present, // Last updated value of sigma_8 (Presently PLANCK+WMAP)
  vnn; //Spectral index of primordial Power spectrum

long N1, N2, N3; // box dimension (grid)
int NF; // Fill every NF grid point 

float LL; // grid spacing in Mpc

long MM; // Number of particles

int zel_flag=1, // memory allocation for zel is 3 times that for nbody
  fourier_flag; //for fourier transfrom
float  DM_m, // Darm matter mass of simulation particle in (10^10 M_sun h^-1) unit 
  norm, // normalize Pk
  pi=M_PI;

io_header header1;

//-----------------------------------------------------------------------------//
//                          from N-body code  done
//-----------------------------------------------------------------------------//



//-----------------------------------------------------------------------------//
//                           needed for N-bdy funcs 
//-----------------------------------------------------------------------------//

float ***ro; // for density
fftwf_plan p_ro; // for FFT
fftwf_plan q_ro; // for FFT

//-----------------------------------------------------------------------------//
//                           needed for N-bdy funcs done  
//-----------------------------------------------------------------------------//

//----------------------arrays for storing ionization data---------------------//

float ***nh, // stores neutral hydrogen on grid points
  ***nhs,    // stores smoothed neutral hydrogen on grid point
  ***ngamma, // stores photon number on grid points
  ***ngammas, // stores smoothed photon number on grid points
  ***nxion;  // stores ionization fractions for different nions on grid points

/*----------------------------GLOBAL VARIABLES DONE----------------------------*/


void main()
{
  long int seed;
  FILE  *inp,*outpp;
  
  int i;
  long ii,jj, kk,ll;
  int sfac;
  
  float vaa;  // final scale factor
  
  //-----------------done variables for non-uniform recombination------------//
  
  float dr,r_min,r_max; // radious for smoothing
  
  char file[100], file1[100], file2[100], num[8], num1[8], num2[8];
  
  float nion,xh1; // to store ionization fraction and neutral hydrogen fraction
  double vion, roion; // to store vol. avg. and mass avg. ionization fraction
  
  int output_flag,in_flag;
  
  long totcluster; // total no. of haloes in halo_catalogue
  float robar,Radii;
  double robarhalo; //no. of dark matter paricle density per (grid)^3
  float vfac; //coefficient for redshift space distortion
  
  float **rra,**vva,**data, //to store particle positions and velocities
    **halo;
  
  double t,T=omp_get_wtime(); // for timing
  
  int Noutput;
  float *nz, *rz, *saveroion ;
  float rl=0.0, rr=0.0; // light speed in km/s
  
  /*---------------------------------------------------------------------------*/
  /* Read input parameters for the simulation from the file "input.nbody_comp" */
  /*---------------------------------------------------------------------------*/
  
  inp=fopen("../../Save/Run/inputs/input.sampling","r");
  fscanf(inp,"%*d%*d");
  fscanf(inp,"%f%f%f%f",&vhh,&vomegam,&vomegalam,&vnn);
  fscanf(inp,"%f%f",&vomegab,&sigma_8_present);
  fscanf(inp,"%d%d%d%d%f",&N1,&N2,&N3,&NF,&LL);
  fscanf(inp,"%*d%*d");
  fscanf(inp,"%*f%*f");  /* time step, final scale factor*/
  fscanf(inp,"%d",&Noutput);
  
  nz=(float*)calloc(Noutput,sizeof(float)); // array to store Noutput 
  rz=(float*)calloc(Noutput,sizeof(float)); // array to store Noutput 
  
  saveroion=(float*)calloc(Noutput,sizeof(float)); 

  for(i=0;i<Noutput;i++)
    fscanf(inp,"%f %f", &nz[i], &rz[i]);
  
  fclose(inp);
  
  //---------------------------------------------------------------------------//
  //-------------parameters read from input file. Check this ------------------//
  
  sfac=8;
  nion=23.21;
  
  //---------------------------------------------------------------------------//
  //system("rm -r ionz_out");
  system("mkdir -p ../../Save/Run/ionz_out_lc");
  
  /*-----------------------------read nbody output-----------------------------*/
  char *filename = "../../Save/Run/inputs/input.lightcone";
  printf(" opening input light cone = %s \n ", filename);
  FILE *fp = fopen(filename, "w");
  fprintf(fp,"%f %f %f %f\n",vhh,vomegam,vomegalam,vomegab);
  printf("%f %f %f %f\n", vhh,vomegam,vomegalam,vomegab);
  fprintf(fp,"%d %d %d\n",N1,N2,N3);
  fprintf(fp,"%f\n",LL);
//Define start
  fprintf(fp,"%d\n",11);
  fprintf(fp,"%d\n",Noutput);

  for(i=0;i<Noutput;i++)
    {
      //-------------------------reading the halo catalogue------------------------//
      
      t=omp_get_wtime();
      strcpy(file1,"../../Save/Run/halo_catalogue_");
      

      sprintf(num1,"%.3f",nz[i]);
      strcat(file1,num1);

      printf(" opening halo_catalogue_  %s\n ", file1);
      
      read_fof(file1,1,&output_flag,&totcluster,halo,&vaa);
      
      halo = allocate_float_2d(totcluster,7);
      
      read_fof(file1,2,&output_flag,&totcluster,halo,&vaa);
      
      printf("ok read halo catalogue = %e\n",omp_get_wtime()-t);
      //printf("%ld\n", MM);
      
      //-------------------------reading the particle data------------------------//
      
      t=omp_get_wtime();
      strcpy(file,"../../Save/Run/sampled_outputs/sampled.nbody_");
      sprintf(num,"%.3f",nz[i]);
      strcat(file,num);
      printf(" opening sampled  %s\n ", file);
      read_sampled(file,1,&seed,&output_flag,&in_flag,rra,vva,&vaa); // only read header
      
      if(i==0)
	{
	  rra = allocate_float_2d(MM,3);
          vva = allocate_float_2d(MM,1);
	  data = allocate_float_2d(MM,5);
	}
      
      read_sampled(file,2,&seed,&output_flag,&in_flag,rra,vva,&vaa); // read data
      
      printf("ok read nbody output = %e\n",omp_get_wtime()-t);
      //printf("%ld\n", MM);
      
      //-----------------------------Redefine grid---------------------------------//
      
      N1=N1/sfac;  N2=N2/sfac;  N3=N3/sfac;// new grid dimensions 
      LL=LL*sfac; 
      robar=MM*8./(1.*N1*N2*N3); // mean numbder density (grid)^{-3}
      
      //---------------------------------------------------------------------------//
      
      for(ii=0;ii<MM;ii++)
	{
	  data[ii][0] = rra[ii][0]/(1.*sfac);
	  data[ii][1] = rra[ii][1]/(1.*sfac);
	  data[ii][2] = rra[ii][2]/(1.*sfac);

	  data[ii][3] = vva[ii][0];
	  
	  data[ii][4] = 8.;  // same mass for all particles
	}

      /*----------------------------------------------------------------*/
      
      for(ii=0;ii<totcluster;ii++)
	{
	  halo[ii][1] /= (1.*sfac);
	  halo[ii][2] /= (1.*sfac);
	  halo[ii][3] /= (1.*sfac);
	}
      
      /*----------------------------------------------------------------*/
      
      if(i==0)
	{
	  Setting_Up_Memory_For_ionz();
	}
      
      /*----------------------------------------------------------------*/
      
      t=omp_get_wtime();
      
      /* calculating the halo mass density at each grid point */

      MM=totcluster;
      cic_vmass(ngamma, halo, 1, 2, 3, 0);  
      
      /*----------------------------------------------------------------*/
      
      MM=header1.npart[1];
      cic_vmass(nh, data, 0, 1, 2, 4);
      
      printf("ok cic_vmass= %e\n",omp_get_wtime()-t);
      
      /*----------------------------------------------------------------*/
      
      t=omp_get_wtime();
      
      for(ii=0;ii<N1;ii++)
	for(jj=0;jj<N2;jj++)
	  for(kk=0;kk<N3;kk++)
	    {
	      if(nh[ii][jj][kk]>nion*ngamma[ii][jj][kk]) // checking ionization condition
		nxion[ii][jj][kk]=nion*ngamma[ii][jj][kk]/nh[ii][jj][kk];  
	      else
		nxion[ii][jj][kk]=1.;
	    }
      
      /*----------------------------------------------------------------*/
      
      //calculating max and min radius for smoothing in units of grid size
      
      r_min=1.;
      r_max=20.0/LL; // Mpc/LL in grid unit
      
      /*----------------------------------------------------------------*/
      /*                        smoothing                               */
      /*----------------------------------------------------------------*/
      
      t=omp_get_wtime();
      
      Radii=r_min;
      
      while(Radii < r_max)
	{
	  for(ii=0;ii<N1;ii++)
	    for(jj=0;jj<N2;jj++)
	      for(kk=0;kk<N3;kk++)
		{
		  nhs[ii][jj][kk]=nh[ii][jj][kk];
		  ngammas[ii][jj][kk]=ngamma[ii][jj][kk];
		}
	  smooth(nhs,Radii);
	  
	  smooth(ngammas,Radii);
	  
	  
	  for(ii=0;ii<N1;ii++)
	    for(jj=0;jj<N2;jj++)
	      for(kk=0;kk<N3;kk++)
		{
		  if(nhs[ii][jj][kk]<=nion*ngammas[ii][jj][kk])  // checking ionization condition
		    nxion[ii][jj][kk]=1.;
		}
	  
	  dr=(Radii*0.1) < 2.0 ? (Radii*0.1) : 2.0; //increment of the smoothing radius
	  Radii += dr;
	}
      
      printf("ok smoothing = %e\n",omp_get_wtime()-t);
      
      /*----------------------------------------------------------------*/
      
      t=omp_get_wtime();
      
      //---------------calculating avg. neutral fraction-------------*/

      sprintf(file2,"%s%.4f","../../Save/Run/ionz_out_lc/HI_map",nz[i]);

      printf(" opening Hi map  %s\n ", file2);
      outpp=fopen(file2,"w");

      fwrite(&N1,sizeof(int),1,outpp);
      fwrite(&N2,sizeof(int),1,outpp);
      fwrite(&N3,sizeof(int),1,outpp);

      roion=0.0;
      
      for(ii=0;ii<N1;ii++)
	for(jj=0;jj<N2;jj++)
	  for(kk=0;kk<N3;kk++)
	    {
  	      xh1=(1.-nxion[ii][jj][kk]);
  	      xh1=(xh1 >0.0)? xh1: 0.0;
	      
  	      nxion[ii][jj][kk]=xh1; // store x_HI instead of x_ion
  	      nhs[ii][jj][kk]=xh1*nh[ii][jj][kk]; // ro_HI on grid
	      
 	      roion+=(double)nhs[ii][jj][kk];
	      fwrite(&nhs[ii][jj][kk],sizeof(float),1,outpp);
            }
	    
      fclose(outpp);
      
      roion/=(1.*N1*N2*N3); // mean HI density
      roion/=robar; // divide by H density to get mass avg. xHI
      
      printf("mass avg. x_HI=%.4f\n", roion);
      printf("%f %f %f\n", nz[i], rz[i], roion);
      fprintf(fp,"%f %f %f\n", nz[i], rz[i], roion);

      
      /*----------------------------------------------------------------*/
      
      density_2_mass(nxion, data, 0, 1, 2, 4);   // get particles HI masses from HI density
      
      sprintf(file2,"%s%.4f%s%.4f","../../Save/Run/ionz_out_lc/HI_part_",nz[i],"_",roion);
      saveroion[i]=roion;
      printf(" opening Hi map part %s\n ", file2);
      outpp=fopen(file2,"w");

      fwrite(&N1,sizeof(int),1,outpp);
      fwrite(&N2,sizeof(int),1,outpp);
      fwrite(&N3,sizeof(int),1,outpp);

      MM=header1.npart[1];
      fwrite(&MM,sizeof(long),1,outpp);
      fwrite(&DM_m,sizeof(float),1,outpp); // DM particle mass in units of 10^10 M_sun/h

      fwrite(&LL,sizeof(float),1,outpp);
      fwrite(&roion,sizeof(double),1,outpp);

      for(ii=0;ii<MM;ii++)
	{
	  fwrite(&data[ii][0],sizeof(float),5,outpp);
	}

      fclose(outpp);
      
      free(halo);
      printf("ok time taken= %e\n\n",omp_get_wtime()-t);
    }

  free(rra);
  free(vva);
    
  printf("done. Total time taken = %dhr %dmin %dsec\n",(int)((omp_get_wtime()-T)/3600), (int)((omp_get_wtime()-T)/60)%60, (int)(omp_get_wtime()-T)%60);


  fclose(fp);
}

