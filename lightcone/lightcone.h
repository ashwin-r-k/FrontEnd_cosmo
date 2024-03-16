float Hf(float aa),Df(float aa),ff(float aa);
float **allocate_float_2d(long N1,int N2);
float  ***allocate_fftwf_3d(long N1,long N2,long N3);
void Setting_Up_Memory_For_ionz();
void cic_vmass(float ***ro_dum, float **data,int xin,int yin,int zin,int min);
void calpow_mom(float ***ro_dum,int Nbin,double* power, double* kmode, double* power_P2,double* power_P4, double* no);
void calpow_mom_k(float ***ro_dum,int Nbin,float kmin,float kmax,double* power, double* kmode, double* power_P2,double* power_P4, double *no);
