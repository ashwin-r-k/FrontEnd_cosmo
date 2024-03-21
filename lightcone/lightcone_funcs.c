#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<fftw3.h>

#include<omp.h>

#include"lightcone.h"


/*---------------------GLOBAL VARIABLES declaration---------------------*/

extern int N1,N2,N3;// box dimension (grid) 
extern int Nbin; // Number of bins to calculate final P(k) (output)

extern long MM, newMM; // Number of particles 
extern float LL; // grid spacing in Mpc

extern float robar;
extern float zc;
extern float pi;

//-------------------arrays for storing data-----------------------------//

extern float ***ro;  // for density or potential
extern fftwf_plan p_ro; // for FFT
extern fftwf_plan q_ro; // for FFT

void Setting_Up_Memory_For_ionz()
{
  //-------------------for multiple  threads---------------------------//
  fftwf_init_threads();
  printf("No of threads = %d\n",omp_get_max_threads()); 
  fftwf_plan_with_nthreads(omp_get_max_threads());
  
  omp_set_num_threads(omp_get_max_threads());
  
  //---------------------done multi thread----------------------------// 
  
  //--------------------arrays for storing data----------------------//

  ro=allocate_fftwf_3d(N1,N2,N3+2);
}

//*************************************************************************

//*************************************************************************

void cic_vmass(float ***ro_dum, float **data, int xin, int yin, int zin, int min)
/* This uses Cloud in Cell for calculating density given posns.*/
/* The i/p is simply the array containing the posns. */
{
  int ii, jj, kk, ix, jy, kz;
  long i, j, k, a[2], b[2], c[2], pin;
  float xx,yy,zz,delx,dely,delz,wx,wy,wz,W;
  
  /* Clear out the array ro. ******/
  for(i=0;i<N1;i++)
    for(j=0;j<N2;j++)
      for(k=0;k<N3;k++)
	ro_dum[i][j][k] = 0.0;
  
  /********************************/
  for(pin=0;pin<MM;pin++)
    { /* begin particle index loop */
      /* (a/b/c)[0] or (a/b/c)[1] can never be greater than (N1/N2/N3) */
      
      a[0]=floor(data[pin][xin]);
      b[0]=floor(data[pin][yin] );
      c[0]=floor(data[pin][zin]);
      
      a[1]=(a[0]+1);
      b[1]=(b[0]+1);
      c[1]=(c[0]+1);
      
      xx=data[pin][xin];
      yy=data[pin][yin];
      zz=data[pin][zin];
      
      /* for each of the 8 corner points */
      for(ii=0;ii<=1;ii++)
	for(jj=0;jj<=1;jj++)
	  for(kk=0;kk<=1;kk++)
	    { /* begin 8 corners loop */
	      delx = xx - a[ii];
	      dely = yy - b[jj];
	      delz = zz - c[kk];
	      
	      ix=a[ii]%N1;
	      jy=b[jj]%N2;
	      kz=c[kk]%N3;
	      
	      /* assigning of weights to the corners */
	      wx=1.0-fabs(delx);
	      wy=1.0-fabs(dely);
	      wz=1.0-fabs(delz);
	      W=wx*wy*wz*data[pin][min]; //multiplying the product of weights with mass of halo
	      
	      ro_dum[ix][jy][kz]+= W;
	    } /* end of <8 grid corners loop>	*/
    } /* end of each particle loop */

} /* end of function cic_ph */


//*************************************************************************

//*************************************************************************

void calpow_mom(float ***ro_dum,int Nbin,double* power, double* kmode, double* power_P2,double* power_P4, double *no)
{ 
  
  /******************** TO FIND POWER SPECTRUM **************************/
  
  long i, j, k, a, b, c, d;
  long index, index1, index2;
  
  double *no2, *no4;
  fftwf_complex *comp_ro;
  
  float fac1, fac2, fac3, m, mu, P2, P4, scale;
  double norml;
  
  norml=1./(1.*N1*N2*N3);
  
  /*************** TAKING FOURIER TRANSFORM OF RO. **************/
  
  fftwf_plan p_ro_dum=fftwf_plan_dft_r2c_3d(N1, N2, N3, &(ro_dum[0][0][0]), (fftwf_complex*)&(ro_dum[0][0][0]),FFTW_ESTIMATE);
  
  fftwf_plan q_ro_dum=fftwf_plan_dft_c2r_3d(N1, N2, N3, (fftwf_complex*)&(ro_dum[0][0][0]), &(ro_dum[0][0][0]),FFTW_ESTIMATE);
  
  
  fftwf_execute(p_ro_dum);
  
  comp_ro=(fftwf_complex*)&(ro_dum[0][0][0]);
  
  /*********** TO FIND POWER SPECTRUM OF RO. **************/
  
  no2=calloc((size_t)Nbin,sizeof(double));
  no4=calloc((size_t)Nbin,sizeof(double));
  
  fac1=1./(1.*N1*N1);
  fac2=1./(1.*N2*N2);
  fac3=1./(1.*N3*N3);
  
  /**************** BINNING POWER SPECTRA **********************/
  
  for(i=0;i<Nbin;++i)
    {
      kmode[i]=0.0;
      
      power[i]=0.0;
      power_P2[i]=0.0;
      power_P4[i]=0.0;
      
      no[i]=0.0;
      no2[i]=0.0;
      no4[i]=0.0;
    }
  
  //***********************************************************************
  
  scale=log10(0.5*N1)/Nbin;
  
  //***********************************************************************

  /*---------------------- BINNING POWER SPECTRA -------------------*/
  
  /*-------------------------- half lines ------------------------- */
  
  for(i=1;i<=N1/2;i++)
    for(j=0;j<=N2/2;j=j+N2/2)
      for(k=0;k<=N3/2;k=k+N3/2)
	{
	  index = i*N2*(N3/2+1) + j*(N3/2+1) + k;
	  
	  m = sqrt(fac1*i*i + fac2*j*j + fac3*k*k);	      
	  
	  mu = (1.0*k)/(N3*m);
	  
	  P2 = 0.5*(3.*mu*mu - 1.);
	  P4 = 0.125*(35.0*powf(mu,4.0)-30.0*mu*mu + 3);
	  
	  d=(int)floorf(log10(m*N1)/scale);  // logarithmic bins
	  
	  if(d>=0 && d<Nbin)
	    {
	      power[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))); 
	      power_P2[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))*P2);
	      power_P4[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))*P4);
	      
	      kmode[d]+=(double)(m);
	      no[d]= no[d] + (double)(1.0);
	      no2[d]=no2[d] + (double)(P2*P2);
	      no4[d]=no4[d] + (double)(P4*P4);
	    }
	}	  
  
  /*----------------------- half planes -----------------------*/
  
  for(i=0;i<N1;i++)
    {
      a=(i>N1/2)? N1-i: i;
      index1 = i*N2*(N3/2+1) ;
      
      for(j=1;j<N2/2;j++) 
	{
	  b=j; 
	  index2 = index1 + j*(N3/2+1) ;
	  
	  for(k=0;k<=N3/2;k=k+N3/2)
	    {
	      c=k;
	      index = index2 + k;
	      
	      m = sqrt(fac1*a*a + fac2*b*b + fac3*c*c);	      
	      
	      mu = (1.0*c)/(N3*m);
	      
	      P2 = 0.5*(3.*mu*mu - 1.);
	      P4 = 0.125*(35.0*powf(mu,4.0)-30.0*mu*mu + 3);
	      
	      d=(int)floorf(log10(m*N1)/scale);//logarithmic bins
	      
	      if(d>=0 && d<Nbin)
		{
		  power[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))); 
		  power_P2[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))*P2);
		  power_P4[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))*P4);
		  
		  kmode[d]+=(double)(m);
		  no[d]= no[d] + (double)(1.0);
		  no2[d]=no2[d] + (double)(P2*P2);
		  no4[d]=no4[d] + (double)(P4*P4);
		}
	      
	    }	  
	}
    }
  
  /**************** half cube **********************/
  
  for(i=0;i<N1;i++)
    {
      a=(i>N1/2)? N1-i: i;
      index1 = i*N2*(N3/2+1) ;
      
      for(j=0;j<N2;j++)
	{
	  b=(j>N2/2)? N2-j: j;
	  index2 = index1 + j*(N3/2+1) ;
	  
	  for(k=1;k<N3/2;k++)
	    {
	      c=k;	  	      
	      index = index2 + k;
	      
	      m = sqrt(fac1*a*a + fac2*b*b + fac3*c*c);	      
	      
	      /* m*(2 * pi/LL) is |k| */
	      /* m=1/2 corresponds to kmode[Nbin-1] i.e. Nyquits */
	      
	      mu = (1.0*c)/(N3*m);
	      
	      P2 = 0.5*(3.*mu*mu - 1.);
	      P4 = 0.125*(35.0*powf(mu,4.0)-30.0*mu*mu + 3);
	      
	      d=(int)floorf(log10(m*N1)/scale);//logarithmic bins
	      
	      if(d>=0 && d<Nbin)
		{
		  power[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))); 
		  power_P2[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))*P2);
		  power_P4[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))*P4);
		  
		  kmode[d]+=(double)(m);
		  no[d]= no[d] + (double)(1.0);
		  no2[d]=no2[d] + (double)(P2*P2);
		  no4[d]=no4[d] + (double)(P4*P4);
		}
	    } /* end k for loop */
	  
	}
    }
  
  //***********************************************************************
  //***********************************************************************
  
  for(i=0;i<Nbin;i++)
    {
      if (no[i]>0.0) 
	{
	  power[i] = pow(LL,3.)*norml*power[i]/(1.0*no[i]);	  
	  kmode[i]=2.*pi*kmode[i]/(no[i]*LL);
	}
      if(no2[i]>0.0)
	power_P2[i] = pow(LL,3.)*norml*power_P2[i]/(1.0*no2[i]);
      
      if(no4[i]>0.0)
	power_P4[i] = pow(LL,3.)*norml*power_P4[i]/(1.0*no4[i]);
    }
  
  //***********************************************************************
  
  for(i=0;i<N1;i++)
    {
      index1 = i*N2*(N3/2+1) ;	  
      for(j=0;j<N2;j++)
	{
	  index2=index1 + j*(N3/2+1) ;
	  for(k=0;k<(N3/2+1);k++)
	    {
	      index=index2 + k;
	      comp_ro[index][0]=comp_ro[index][0]*norml;
	      comp_ro[index][1]=comp_ro[index][1]*norml;
	    }
	}
    }
  
  /*  now convert the array back to real space */
  
  //***********************************************************************
  
  fftwf_execute(q_ro_dum);
  
  fftwf_destroy_plan(p_ro_dum);
  fftwf_destroy_plan(q_ro_dum); 
  
} /* end function */

//*************************************************************************



//*************************************************************************

void calpow_mom_k(float ***ro_dum,int Nbin,float kmin,float kmax,double* power, double* kmode, double* power_P2,double* power_P4, double *no)
{ 
  
  /******************** TO FIND POWER SPECTRUM **************************/
  
  long i, j, k, a, b, c, d;
  long index, index1, index2;
  
  double *no2, *no4;
  fftwf_complex *comp_ro;
  
  float fac1, fac2, fac3, m, mu, P2, P4, scale;
  double norml;
  
  float tpibyL;
  tpibyL = 2.0*pi/LL; // 2 pi /LL
  
  norml=1./(1.*N1*N2*N3);
  
  /*************** TAKING FOURIER TRANSFORM OF RO. **************/
  
  fftwf_plan p_ro_dum=fftwf_plan_dft_r2c_3d(N1, N2, N3, &(ro_dum[0][0][0]), (fftwf_complex*)&(ro_dum[0][0][0]),FFTW_ESTIMATE);
  
  fftwf_plan q_ro_dum=fftwf_plan_dft_c2r_3d(N1, N2, N3, (fftwf_complex*)&(ro_dum[0][0][0]), &(ro_dum[0][0][0]), FFTW_ESTIMATE);
  
  
  fftwf_execute(p_ro_dum);
  
  comp_ro=(fftwf_complex*)&(ro_dum[0][0][0]);
  
  /*********** TO FIND POWER SPECTRUM OF RO. **************/
  
  no2=calloc((size_t)Nbin,sizeof(double));
  no4=calloc((size_t)Nbin,sizeof(double));
  
  fac1=1./(1.*N1*N1);
  fac2=1./(1.*N2*N2);
  fac3=1./(1.*N3*N3);
  
  /**************** BINNING POWER SPECTRA **********************/
  
  for(i=0;i<Nbin;++i)
    {
      kmode[i]=0.0;
      
      power[i]=0.0;
      power_P2[i]=0.0;
      power_P4[i]=0.0;
      
      no[i]=0.0;
      no2[i]=0.0;
      no4[i]=0.0;
    }
  
  //***********************************************************************
  
  scale=(log10(kmax)-log10(kmin))/Nbin;
  
  //***********************************************************************
  
  /************************ BINNING POWER SPECTRA *************************/
  
  /*-------------------------- half lines ------------------------- */
  
  for(i=1;i<=N1/2;i++)
    for(j=0;j<=N2/2;j=j+N2/2)
      for(k=0;k<=N3/2;k=k+N3/2)
	{
	  a=i;
	  b=j;
	  c=k;
	  
	  index = i*N2*(N3/2+1) + j*(N3/2+1) + k;
	  
	  m = tpibyL*sqrt(fac1*a*a + fac2*b*b + fac3*c*c);  /* m is |k| */      
	  
	  mu = (tpibyL*c)/(N3*m);
	  
	  P2 = 0.5*(3.*mu*mu- 1.);
	  P4 = 0.125*(35.0*powf(mu,4.0)-30.0*mu*mu + 3);
	  
	  d=(int)floorf((log10(m)-log10(kmin))/scale);  //logarithmic bins
	  
	  if(d>=0 && d<Nbin)
	    {
	      power[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))); 
	      power_P2[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))*P2);
	      power_P4[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))*P4);
	      
	      kmode[d]+=(double)(m);
	      no[d]= no[d] + (double)(1.0);
	      no2[d]=no2[d] + (double)(P2*P2);
	      no4[d]=no4[d] + (double)(P4*P4);
	    }
	}	  
  
  /*----------------------- half planes -----------------------*/
  
  for(i=0;i<N1;i++)
    {
      a=(i>N1/2)? N1-i: i;
      index1 = i*N2*(N3/2+1) ;
      
      for(j=1;j<N2/2;j++) 
	{
	  b=j; 
	  index2 = index1 + j*(N3/2+1) ;
	  
	  for(k=0;k<=N3/2;k=k+N3/2)
	    {
	      c=k;
	      index = index2 + k;
	      
	      m = tpibyL*sqrt(fac1*a*a + fac2*b*b + fac3*c*c);  /* m is |k| */      
	      
	      mu = (tpibyL*c)/(N3*m);
	      
	      P2 = 0.5*(3.*mu*mu- 1.);
	      P4 = 0.125*(35.0*powf(mu,4.0)-30.0*mu*mu + 3);
	      
	      d=(int)floorf((log10(m)-log10(kmin))/scale);  //logarithmic bins
	      
	      if(d>=0 && d<Nbin)
		{
		  power[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))); 
		  power_P2[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))*P2);
		  power_P4[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))*P4);
		  
		  kmode[d]+=(double)(m);
		  no[d]= no[d] + (double)(1.0);
		  no2[d]=no2[d] + (double)(P2*P2);
		  no4[d]=no4[d] + (double)(P4*P4);
		}
	      
	    }	  
	}
    }
  
  /**************** half cube **********************/
  
  for(i=0;i<N1;i++)
    {
      a=(i>N1/2)? N1-i: i;
      index1 = i*N2*(N3/2+1) ;
      
      for(j=0;j<N2;j++)
	{
	  b=(j>N2/2)? N2-j: j;
	  index2 = index1 + j*(N3/2+1) ;
	  
	  for(k=1;k<N3/2;k++)
	    {
	      c=k;	  	      
	      index = index2 + k;
	      
	      m = tpibyL*sqrt(fac1*a*a + fac2*b*b + fac3*c*c);  /* m is |k| */      
	      
	      mu = (tpibyL*c)/(N3*m);
	      
	      P2 = 0.5*(3.*mu*mu- 1.);
	      P4 = 0.125*(35.0*powf(mu,4.0)-30.0*mu*mu + 3);
	      
	      d=(int)floorf((log10(m)-log10(kmin))/scale);  //logarithmic bins
	      
	      if(d>=0 && d<Nbin)
		{
		  power[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))); 
		  power_P2[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))*P2);
		  power_P4[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))*P4);
		  
		  kmode[d]+=(double)(m);
		  no[d]= no[d] + (double)(1.0);
		  no2[d]=no2[d] + (double)(P2*P2);
		  no4[d]=no4[d] + (double)(P4*P4);
		}
	    } /* end k for loop */
	  
	}
    }
  
  //***********************************************************************
  //***********************************************************************
  
  for(i=0;i<Nbin;i++)
    {
      if (no[i]>0.0) 
	{
	  power[i] = pow(LL,3.)*norml*power[i]/(1.0*no[i]);	  
	  kmode[i] = kmode[i]/(no[i]*1.0);
	}
      if(no2[i]>0.0)
	power_P2[i] = pow(LL,3.)*norml*power_P2[i]/(1.0*no2[i]);
      
      if(no4[i]>0.0)
	power_P4[i] = pow(LL,3.)*norml*power_P4[i]/(1.0*no4[i]);
    }
  
  //***********************************************************************
  
  for(i=0;i<N1;i++)
    {
      index1 = i*N2*(N3/2+1) ;	  
      for(j=0;j<N2;j++)
	{
	  index2=index1 + j*(N3/2+1) ;
	  for(k=0;k<(N3/2+1);k++)
	    {
	      index=index2 + k;
	      comp_ro[index][0]=comp_ro[index][0]*norml;
	      comp_ro[index][1]=comp_ro[index][1]*norml;
	    }
	}
    }
  
  /*  now convert the array back to real space */
  
  //***********************************************************************
  
  fftwf_execute(q_ro_dum);
  
  fftwf_destroy_plan(p_ro_dum);
  fftwf_destroy_plan(q_ro_dum); 
  
} /* end function */
