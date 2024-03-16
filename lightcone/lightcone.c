#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<fftw3.h>
#include<omp.h>

#include"lightcone.h"

int N1, N2, N3; // box dimension (grid) 
int Nbin=10, sfac=8; // Number of bins to calculate final P(k) (output) // check

float LL; // grid spacing in Mpc
float DM_m;
float vhh,vomegam,vomegalam,vomegab;

long MM, newMM=0; // Number of particles

float robar; // check
float zc; // check

float pi=M_PI;

float ***ro; // for density
fftwf_plan p_ro; // for FFT
fftwf_plan q_ro; // for FFT

float *rz, *nionz;

void main()
{
  FILE  *inp,*outpp;
  long ii, jj, kk, ll, index;
  
  char file[100], num[8];
  
  float **data, dummy; //to store particle positions and velocities
  
  double t,T=omp_get_wtime(); // for timing
  
  int Noutput, Ninput;
  int start; // check
  float *nz, inR, finR;
  
  /*---------------------------------------------------------------------------*/
  /* Read input parameters for the simulation from the file "input.nbody_comp" */
  /*---------------------------------------------------------------------------*/
  
  inp=fopen("../../Save/Run/inputs/input.lightcone","r");

  fscanf(inp,"%f %f %f %f",&vhh,&vomegam,&vomegalam,&vomegab);
  fscanf(inp,"%d %d %d",&N1,&N2,&N3);
  fscanf(inp,"%f",&LL);
  printf("N=%d, LL=%f\n",N1,LL);
  
  fscanf(inp,"%d",&start);
  fscanf(inp,"%d",&Noutput);

  nz=(float*)calloc(Noutput,sizeof(float)); // array to store Noutput 
  rz=(float*)calloc(Noutput,sizeof(float));  
  nionz=(float*)calloc(Noutput,sizeof(float));
  
  for(ii=0;ii<Noutput;ii++)
    fscanf(inp,"%f %f %f",&nz[ii],&rz[ii],&nionz[ii]);
  
  fclose(inp);

  /*---------------------------------------------------------------------------*/
  
  double *power_P0, *power_P2, *power_P4, *kmode; // arrays for power spectrum 
  double *no;
  
  Setting_Up_Memory_For_ionz();
  
  kmode=calloc((size_t)Nbin,sizeof(double));
  power_P0=calloc((size_t)Nbin,sizeof(double));
  power_P2=calloc((size_t)Nbin,sizeof(double));
  power_P4=calloc((size_t)Nbin,sizeof(double)); 
  no=calloc((size_t)Nbin,sizeof(double));
  
  double tmpRo, *roion, *vion; // to store mass avg. ionization fraction and Mean Tb
  
  roion=calloc((size_t)(N3),sizeof(double));
  vion =calloc((size_t)(N3),sizeof(double));
  /*---------------------------------------------------------------------------*/
  
  inR = rz[start]; // check
  finR= inR + LL*N3; // check

  printf("inR=%f, finR=%f in Mpc\n", inR, finR);

  float theta_max=LL*N1*1./(2.*finR);
  float delta_theta=LL/finR;
  printf("theta_max=%f, delta_theta=%f\n", theta_max, delta_theta);

  float vfac, c = 299792.458; // light speed in km/sec
  float nu_i, nu_f, nu, nu_e = 1420.405751;
  nu_i = nu_e/(1.+nz[start]);
  printf("nui=%f in MHz.\n", nu_i);
  
  /*---------------------------------------------------------------------------*/
  
  Ninput=start;
  while(rz[Ninput]<finR)
    Ninput++;
  printf("%d no. of outputs are needed to create the lightcone.\n", Ninput - start);
  
  outpp=fopen("lc_part","w");
  
  /*---------------------------------------------------------------------------*/

  
  double r1, r2, RR, tmpR;
  r2 = rz[start];
    
  for(ll=start;ll<Ninput;ll++)
    {
      sprintf(file,"%s%.4f%s%.4f","../../Save/Run/ionz_out_lc/HI_part_",nz[ll],"_",nionz[ll]);
      printf("ll = %d  Ninp = %d\n",ll,Ninput );
      printf("File trying %s\n",file );
      inp=fopen(file,"r");

      printf("File opened %s\n",file );

      /*-----------------------------------------------------------------------*/
      fread(&N1,sizeof(int),1,inp);
      fread(&N2,sizeof(int),1,inp);
      fread(&N3,sizeof(int),1,inp);

      fread(&MM,sizeof(long),1,inp);
      fread(&DM_m,sizeof(float),1,inp);

      fread(&LL,sizeof(float),1,inp);
      fread(&tmpRo,sizeof(double),1,inp);
      
      if(ll==start)
	data = allocate_float_2d(MM,5);
      
      for(ii=0;ii<MM;ii++)
	{
	  fread(&data[ii][0],sizeof(float),5,inp);
	  
	  data[ii][0] = data[ii][0]*LL - .5*N1*LL; // in Mpc 
	  data[ii][1] = data[ii][1]*LL - .5*N2*LL; // in Mpc 
	  data[ii][2] = data[ii][2]*LL + inR; // in Mpc
	  
	  RR = sqrt(data[ii][0]*data[ii][0] + data[ii][1]*data[ii][1] + data[ii][2]*data[ii][2]);

	  data[ii][0] = data[ii][0]/RR;
	  data[ii][1] = data[ii][1]/RR;
	  data[ii][2] = RR;
	}
      fclose(inp);

      /*-----------------------------------------------------------------------*/
      
      vfac = Hf(1./(1. + nz[ll]))*100.*vhh/(c*(1. + nz[ll])); // H in (km/sec)/Mpc and c in km/sec

      r1=r2;
      if(rz[ll+1]<finR)
	r2 = rz[ll+1];
      else
	{
	  r2 = finR;
	  nu_f = nu_e*(1. - vfac*(r2 - r1))/(1. + nz[ll]);
	  printf("nu_f=%f in MHz. last step done.\n", nu_f);
	}      
      

      for(ii=0;ii<MM;ii++)
        {
  	  if(data[ii][2]>=r1 && data[ii][2]<r2)
  	    {
	      if(fabs(data[ii][0])<=theta_max && fabs(data[ii][1])<=theta_max)
		{
		  data[ii][3] = nu_e*(1. - vfac*(data[ii][2] - rz[ll]) - data[ii][3])/(1. + nz[ll]);
		  
		  fwrite(&data[ii][0],sizeof(dummy),5,outpp);
		  newMM++;
		}
  	    }
  	}
      
    }
  fclose(outpp);
  free(data);
  
  /*---------------------------------------------------------------------------*/
  
  printf("MM=%ld, newMM=%ld\n",MM, newMM);

  data = allocate_float_2d(newMM,5);

  float nu_cuti=0.0, nu_cutf=0.0;
  
  inp=fopen("lc_part","r");
  
  for(ii=0;ii<newMM;ii++)
    {
      fread(&data[ii][0],sizeof(float),5,inp);

      if((data[ii][3] - nu_i) > nu_cuti)
	nu_cuti = data[ii][3] - nu_i; 
	
      if((nu_f - data[ii][3]) > nu_cutf)
	nu_cutf = nu_f - data[ii][3];
    }
  fclose(inp);
  printf("nu_cuti=%f, nu_cutf=%f in MHz.\n", nu_cuti, nu_cutf);
    
  /*---------------------------------------------------------------------------*/
  /*-------------------------- check before RSD -------------------------------*/

  /* for(ii=0;ii<newMM;ii++) */
  /*   { */
  /*     data[ii][0] = (data[ii][0] + theta_max)/delta_theta; // in Mpc */
  /*     data[ii][1] = (data[ii][1] + theta_max)/delta_theta; // in Mpc */
      
  /*     data[ii][2] = (data[ii][2] - inR)/LL; */
  /*   } */

  /* robar=MM/(1.*N1*N2*N3); */
  /* printf("%f\n", robar); */

  /* MM=newMM; */
  /* printf("%ld\n", MM); */
  /* cic_vmass(ro, data, 0, 1, 2, 4);  // convert particles HI masses to  HI density   */
  
  /* for(ii=0;ii<N1;ii++) */
  /*   for(jj=0;jj<N2;jj++) */
  /*     for(kk=0;kk<N3;kk++) */
  /* 	{ */
  /* 	  vion[kk]+=(double)ro[ii][jj][kk]; */
  /* 	}   */
  
  /* for(kk=0;kk<N3;kk++) */
  /*   { */
  /*     vion[kk] = vion[kk]/(1.*N1*N2); // mean HI density */
      
  /*     vion[kk]/=robar; // divide by H density to get mass avg. xHI */
      
  /*     printf("mass avg. x_HI = %f\n", vion[kk]); */
  /*   } */

  /*---------------------------------------------------------------------------*/
  /*---------------------------------------------------------------------------*/
  /*---------------------------- chop -----------------------------------------*/

  float old_delta_nu;
  old_delta_nu = (nu_i - nu_f)/(1.*N3);
  printf("Bandwidth=%f in MHz. delta_nu=%f\n", nu_i - nu_f, old_delta_nu);

  float delta_nu, nu_c;
  nu_i = nu_i - nu_cuti;
  nu_f = nu_f + nu_cutf;
  
  delta_nu = (nu_i - nu_f)/(1.*N3);
  nu_c = nu_i - .5*delta_nu*N3;
  zc = nu_e/nu_c - 1.;
  printf("After RSD, bandwidth is reduced to %f in MHz at nu_c=%f or zc=%f.\nThe delta_nu=%f in MHz\n", nu_i - nu_f, nu_c, zc, delta_nu);
  printf("After RSD, z_i=%f and z_f=%f.\n", nu_e/nu_i - 1., nu_e/nu_f - 1.);
  
  /*------------------------------ check --------------------------------------*/
  float mass_fac=1., rc;  // question
  mass_fac = LL*LL*LL;
  //mass_fac = 0.508789/LL;  // 0.508789 is the delta_r after RSD
  mass_fac = mass_fac/(delta_theta*delta_theta*delta_nu);
  mass_fac = mass_fac*vhh*100.*nu_e/(c);
  
  rc = (finR + inR)*.5;
  //printf("rc=%f in Mpc.\n", rc);
  printf("mass fraction at r_f = %f.\n", mass_fac/(finR*finR));
  /*---------------------------------------------------------------------------*/

  MM=0;
  outpp=fopen("lc_part_rs","w");

  for(ii=0;ii<newMM;ii++)
    {
      if(data[ii][3]<=nu_i && data[ii][3]>=nu_f)
	{
	  data[ii][4] = data[ii][4]*mass_fac/(data[ii][2]*data[ii][2]); // check
	  //data[ii][4] = data[ii][4]*mass_fac/(finR*finR); // check
	  
	  data[ii][0] = (data[ii][0] + theta_max)/delta_theta;
	  data[ii][1] = (data[ii][1] + theta_max)/delta_theta;
	  data[ii][3] = (nu_i - data[ii][3])/delta_nu;
	  
	  fwrite(&data[ii][0],sizeof(dummy),5,outpp);
	  MM++;
	}
    }
  fclose(outpp);
  free(data);

  /*---------------------------------------------------------------------------*/

  //printf("%d %ld %ld %f %f\n", 8*N1*N2*N3, newMM, MM, (8.*N1*N2*N3 - newMM)/(8.*N1*N2*N3), (8.*N1*N2*N3 - 1.*MM)/(8.*N1*N2*N3)); // question
  //robar = 8.;
  //robar=newMM/(1.*N1*N2*N3);
  robar=MM/(1.*N1*N2*N3);
  printf("robar=%f\n", robar); // question

  data = allocate_float_2d(MM,5);

  inp=fopen("lc_part_rs","r");
  
  for(ii=0;ii<MM;ii++)
    fread(&data[ii][0],sizeof(float),5,inp);
  
  fclose(inp);

  /*---------------------------------------------------------------------------*/
  
  cic_vmass(ro, data, 0, 1, 3, 4);  // convert particles HI masses to  HI density

  sprintf(file,"%s%.4f","lc_maprs_",zc);
  outpp=fopen(file,"w");
      
  fwrite(&N1,sizeof(int),1,outpp);
  fwrite(&N2,sizeof(int),1,outpp);
  fwrite(&N3,sizeof(int),1,outpp);
            
  for(ii=0;ii<N1;ii++)
    for(jj=0;jj<N2;jj++)
      for(kk=0;kk<N3;kk++)
	{
	  roion[kk]+=1.*ro[ii][jj][kk];
	  ro[ii][jj][kk] = 140.*ro[ii][jj][kk]*vomegab*vhh/robar;
	  fwrite(&ro[ii][jj][kk],sizeof(float),1,outpp);
	  
	  vion[kk]+=1.*ro[ii][jj][kk];
	}
  
  fclose(outpp);
  
  /*----------------------------------------------------------------*/
  
  calpow_mom(ro,Nbin,power_P0,kmode,power_P2,power_P4,no); // calculates moments of redshift space power spectrum
  
  sprintf(file,"%s%.4f","pk.lcrs_",zc);
  outpp=fopen(file,"w");
  
  for(ii=0;ii<Nbin;ii++)
    fprintf(outpp,"%e %e %e %e %ld\n",kmode[ii],power_P0[ii],power_P2[ii],power_P4[ii],(long)no[ii]);
  
  fclose(outpp);
  
  /*----------------------------------------------------------------*/
  sprintf(file,"%s%.4f","z_xHI_",zc);
  outpp=fopen(file,"w");
  
  for(kk=0;kk<N3;kk++)
    {
      roion[kk] = roion[kk]/(1.*N1*N2); // mean HI density
      vion[kk]  = vion[kk]/(1.*N1*N2); // mean Tb in mK
      
      roion[kk] = roion[kk]*(inR + kk*LL)*(inR + kk*LL)/(mass_fac);
      //roion[kk] = roion[kk]*finR*finR/(mass_fac);
      
      roion[kk]/=robar; // divide by H density to get mass avg. xHI
      fprintf(outpp, "%f %f\n", nu_e/(nu_i - kk*delta_nu) - 1., roion[kk]);
      printf("avg. Tb=%f mk\n",vion[kk]);
    }
  fclose(outpp);
  /*----------------------------------------------------------------*/

  printf("done. Total time taken = %dhr %dmin %dsec\n",(int)((omp_get_wtime()-T)/3600), (int)((omp_get_wtime()-T)/60)%60, (int)(omp_get_wtime()-T)%60);
}
