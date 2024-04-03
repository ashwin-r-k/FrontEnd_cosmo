def Zarray(start,stop,step):
    list_Z=[]
    i=start
    while i <= stop :
        list_Z.append(round(i,5))
        i=i+step
    return list_Z
    
seed = -100012
Nbin = 10
hh = 0.6704
omega_m = 0.3183
omega_l = 0.6817
spectral_index =0.9619
omega_baryon = 0.04902
sigma_8 = 0.8347
box=pow(2,7)
fraction_fill = 2
LL= 0.14
a_initial = 0.008
delta_a= 0.004
output_flag=0
pk_flag=1
redshift_values = Zarray(5,10,1)[::-1]
Len_redshift = len(redshift_values)


