#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>

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
int NF, // Fill every NF grid point 
  Nbin; // Number of bins to calculate final P(k) (output)

float LL; // grid spacing in Mpc

long MM; // Number of particles

int zel_flag=1, // memory allocation for zel is 3 times that for nbody
  fourier_flag; //for fourier transfrom
float  DM_m, // Darm matter mass of simulation particle in (10^10 M_sun h^-1) unit 
  norm, // normalize Pk
  pi=M_PI;

float *rz;

io_header header1;

int write_sampled(char *fname, long int seed, int output_flag, float **rra, float **vva, float vaa, float rc);
float **allocate_float_2d(long N1,int N2);
int read_output(char *fname, int read_flag,long int *seed,int *output_flag,int *in_flag,float **rra,float **vva,float *aa);

//-----------------------------------------------------------------------------//
//                          from N-body code  done
//-----------------------------------------------------------------------------//


void main()
{
  long int seed;
  FILE  *inp,*outpp;
  
  int i;
  long ii,jj, kk,ll, tmp;
  int sfac;
  
  float vaa;  // final scale factor
  
  //-----------------done variables for non-uniform recombination------------//
  
  double *power_P0, *power_P2, *power_P4, *kmode; // arrays for power spectrum 
  double *no;
  
  float aa; 
  
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
    
  int Noutput, oflag;
  float *nz, rc;
  
  long *index;
  long new_MM;
  
  /*---------------------------------------------------------------------------*/
  /* Read input parameters for the simulation from the file "input.nbody_comp" */
  /*---------------------------------------------------------------------------*/
  
  inp=fopen("../../Save/Run/inputs/input.sampling","r");
  fscanf(inp,"%ld%*d",&tmp);
  fscanf(inp,"%*f%*f%*f%*f");
  fscanf(inp,"%*f%*f");
  fscanf(inp,"%ld%ld%ld%*d%*f",&tmp,&tmp,&tmp);
  fscanf(inp,"%d%*d", &oflag);
  fscanf(inp,"%*f%*f");  /* time step, final scale factor*/
  fscanf(inp,"%d",&Noutput);

  nz=(float*)calloc(Noutput,sizeof(float)); // array to store Noutput
  rz=(float*)calloc(Noutput,sizeof(float)); // array to store Noutput
  
  for(i=0;i<Noutput;i++)
    fscanf(inp,"%f %f",&nz[i],&rz[i]);
  
  fclose(inp);
  
  /*-----------------------------read nbody output-----------------------------*/
  system("mkdir -p ../../Save/Run/sampled_outputs");
  
  for(i=0;i<Noutput;i++)
    {
      strcpy(file,"../../Save/Run/output.nbody_");
      sprintf(num,"%.3f",nz[i]);
      strcat(file,num);

      printf("%s\n", file);

      read_output(file,1,&seed,&output_flag,&in_flag,rra,vva,&vaa); // only read header
      if(i==0)
	{
	  rra = allocate_float_2d(MM,3);
          vva = allocate_float_2d(MM,3);
	}
      
      read_output(file,2,&seed,&output_flag,&in_flag,rra,vva,&vaa); // read data
      
      aa = vaa;
      rc = rz[i];
      MM = MM/8;
      DM_m = DM_m*8.;
      
      printf("%ld %ld\n", MM*8, MM);
      
      strcpy(file,"../../Save/Run/sampled_outputs/sampled.nbody_");
      sprintf(num,"%.3f",nz[i]);
      strcat(file,num);
      
      write_sampled(file, seed, oflag, rra, vva, aa, rc);

      /* for(ii=0; ii<8; ii++) */
      /* 	{ */
      /* 	  srand(12345+ii); */
      /* 	  index[ii] = ii*8 + (long)round(8.*rand()/(RAND_MAX+1.0)); */
      /* 	  printf("%ld\n", index[ii]); */
      /* 	} */
    }


}




/*-------------------------------------------------------------------------------------------------------*/

int write_sampled(char *fname, long int seed, int output_flag, float **rra, float **vva, float vaa, float rc)
{
  FILE *fp1; 
  long   ii, jj; 
  fp1=fopen(fname,"w");   
  int dummy, kk; 
  // set header structure values 
  
  for(kk=0;kk<6;++kk) 
    { 
      header1.npart[kk]=0; 
      header1.mass[kk]=0.; 
      header1.npartTotal[kk]=0; 
    } 
  header1.npart[1]=(long)MM; //DM particles in this snapshot file 
  header1.mass[1]=DM_m; //DM particle mass in units of 10^10 M_sun/h 
  header1.npartTotal[1]=(long)MM; //Total DM particles in simulation 
  header1.time=vaa;     //scale factor of  nbody output     
  header1.redshift=1./vaa-1.; 
  header1.flag_sfr=0; 
  header1.flag_cooling=0; 
  header1.flag_feedback=0; 
  header1.num_files=1;  //no. of files in each snapshot 
  header1.BoxSize=N1*LL*1000*vhh;//simulation box size in kpc/h 
  header1.Omega0=vomegam;    //Omega_m   
  header1.OmegaLambda=vomegalam;     //OmegaLambda     
  header1.HubbleParam=vhh;     //HubbleParam 
  header1.Omegab=vomegab;    //Omega_b      
  header1.sigma_8_present=sigma_8_present;  
  header1.Nx=N1;
  header1.Ny=N2;
  header1.Nz=N3; 
  header1.LL=LL; 
  header1.output_flag=output_flag; 
  header1.in_flag=zel_flag; 
  header1.seed=seed; 
  // done setting header  
  if(output_flag!=1) 
    { 
      //final paticle positions stored in kpc/h unit and velocity written in km/sec unit  
      for (ii=0;ii<MM*8;++ii) 
 	{	   
 	  rra[ii][0]=rra[ii][0]*LL*1000.*vhh;//coordinates in kpc/h 
 	  rra[ii][1]=rra[ii][1]*LL*1000.*vhh; 
 	  rra[ii][2]=rra[ii][2]*LL*1000.*vhh; 
	  
 	  vva[ii][0]=vva[ii][0]*LL*vhh*100./vaa ; //peculiar velocities in km/sec 
 	  vva[ii][1]=vva[ii][1]*LL*vhh*100./vaa; 
 	  vva[ii][2]=vva[ii][2]*LL*vhh*100./vaa; 
 	}
      printf("ok unit change\n");
    }
  
  // writing header
  
  fwrite(&dummy,sizeof(dummy),1,fp1);     
  fwrite(&header1,sizeof(io_header),1,fp1);     
  fwrite(&dummy,sizeof(dummy),1,fp1); 

  float vfac, tmp1;
  double RR, tmpR, VV;
  float cv = 299792.458; // light speed in km/sec
  float nu; 
  
  rc = rc*1000.*vhh; // in kpc/h unit
  //printf("%f\n", rc);
  
  // writing data  
  fwrite(&dummy,sizeof(dummy),1,fp1); 
  for(ii=0;ii<MM;ii++)
    {
      srand(12345+ii);
      jj = ii*8 + (long)round(8.*rand()/(RAND_MAX+1.0));
      
      fwrite(&rra[jj][0],sizeof(float),3,fp1);

      RR = (rra[jj][0] - .5*N1*LL*1000.*vhh)*(rra[jj][0] - .5*N1*LL*1000.*vhh) + (rra[jj][1] - .5*N2*LL*1000.*vhh)*(rra[jj][1] - .5*N2*LL*1000.*vhh);
      tmpR = (rra[jj][2] + rc);
      RR = sqrt(RR + tmpR*tmpR);

      VV = (rra[jj][0] - .5*N1*LL*1000.*vhh)*vva[jj][0] + (rra[jj][1] - .5*N2*LL*1000.*vhh)*vva[jj][1];
      VV = (VV + tmpR*vva[jj][2])/RR;

      nu = VV/cv;

      //nu = nu_e*vaa*(1. - vfac*(RR - rc)/cv - VV/cv);
      //tmp1 = rra[jj][2]/(LL*1000.*vhh); 
      //if(tmp1<1. || tmp1>511.)
      //printf("%f %lf %f\n", tmp1, nu, nu);
      
      fwrite(&nu,sizeof(float),1,fp1);
    }
  
  fwrite(&dummy,sizeof(dummy),1,fp1); 

  fclose(fp1);


  if(output_flag!=1) 
    { 
      //paticle positions restored in grid units  
      for (ii=0;ii<MM*8;++ii) 
 	{	   
 	  rra[ii][0]=rra[ii][0]/(LL*1000.*vhh);  //coordinates in kpc/h 
 	  rra[ii][1]=rra[ii][1]/(LL*1000.*vhh); 
 	  rra[ii][2]=rra[ii][2]/(LL*1000.*vhh); 
	  
 	  vva[ii][0]=(vva[ii][0]*vaa)/(LL*vhh*100.);  //peculiar velocities in km/sec 
 	  vva[ii][1]=(vva[ii][1]*vaa)/(LL*vhh*100.);
 	  vva[ii][2]=(vva[ii][2]*vaa)/(LL*vhh*100.);
 	} 
    }

} 

float **allocate_float_2d(long N1,int N2) /* Typically N1 is number of particles and N2 is number of dimensions namely 3 */
{
  float **xxa, *xx;
  long ii;

  xxa=(float**)malloc(N1 *  sizeof(float*));
  if(!(xx = (float *) calloc((size_t)(N1*N2),sizeof(float))))
    {
      printf("error in allocate_float_2d");
      exit(0);
    }

  for(ii=0;ii<N1;++ii)
    xxa[ii]=xx + N2*ii ;

return(xxa);
}


int read_output(char *fname, int read_flag,long int *seed,int *output_flag,int *in_flag,float **rra,float **vva,float *aa)
{
  FILE *fp1;
  long ii;
  float vaa;
  fp1=fopen(fname,"r");
  
  int kk,dummy;
  
  // header reading
  fread(&dummy,sizeof(dummy),1,fp1);
  fread(&header1,sizeof(io_header),1,fp1);
  fread(&dummy,sizeof(dummy),1,fp1);
  
  
  vaa=(float)header1.time;     //scale factor of  nbody output
  *aa=vaa;
  
  MM=(long)header1.npartTotal[1]; //Total DM particles in this simulation
  DM_m=(float)header1.mass[1]; //DM particle mass in units of 10^10 M_sun/h
  
  vomegam=(float)header1.Omega0;    //Omega_
  vomegalam=(float)header1.OmegaLambda;     //OmegaLambda
  vhh=(float)header1.HubbleParam;     //HubbleParam
  vomegab=(float)header1.Omegab;    //Omega_b
  sigma_8_present=(float)header1.sigma_8_present;
  N1=(long)header1.Nx;
  N2=(long)header1.Ny;
  N3=(long)header1.Nz;
  LL=(float)header1.LL;
  *output_flag=header1.output_flag; // input units ? (!=1) => kp/h; else grid
  *in_flag=header1.in_flag; //  input file generated by ? 1 => zel  else nbody
  *seed=header1.seed;
  
  if(read_flag!=1)
    {
      fread(&dummy,sizeof(dummy),1,fp1);
      for(ii=0;ii<MM;ii++)
	{
	  fread(&rra[ii][0],sizeof(float),3,fp1);
	  
	  if(*output_flag!=1)
	    {
	      rra[ii][0]=rra[ii][0]/(LL*1000.*vhh); //coordinates in grid
	      rra[ii][1]=rra[ii][1]/(LL*1000.*vhh);
	      rra[ii][2]=rra[ii][2]/(LL*1000.*vhh);
	      
	      rra[ii][0] = rra[ii][0]-1.0*N1*(long)(floor(rra[ii][0])/(1.*N1));
	      rra[ii][1] = rra[ii][1]-1.0*N2*(long)(floor(rra[ii][1])/(1.*N2));  // imposing periodic boundary condition
	      rra[ii][2] = rra[ii][2]-1.0*N3*(long)(floor(rra[ii][2])/(1.*N3));
	    }
	}
      fread(&dummy,sizeof(dummy),1,fp1);
      
      fread(&dummy,sizeof(dummy),1,fp1);
      for(ii=0;ii<MM;ii++)
	{
	  fread(&vva[ii][0],sizeof(float),3,fp1);
	  
	  if (*output_flag!=1)
	    {
	      vva[ii][0]=vva[ii][0]/(LL*vhh*100./vaa); //velocities in 
	      vva[ii][1]=vva[ii][1]/(LL*vhh*100./vaa);
	      vva[ii][2]=vva[ii][2]/(LL*vhh*100./vaa);
	    }
	}
      fread(&dummy,sizeof(dummy),1,fp1);
    }
  
  fclose(fp1);
}
	  
//*************************************************************************
