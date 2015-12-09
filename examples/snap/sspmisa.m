 iopt=[  2  1  0  0  0  0  1  0  0  0  0  0  1  1  0  1  0  0  1  1  0  0  0  0  0  0  0  0  0  1];
 freq=[
       250.0000
 ];
 nplot_c= 1;
 npx( 1)=   51;
 npy( 1)=   51;
 yp( 1,:)='linspace(    1477.500000,    1482.500000,  51)';
 xp( 1,:)='linspace(    1497.500000,    1502.500000,  51)';
 ncx( 1)=   23;
title_x( 1,:)='Water sound speed (m/s)                 ';
 ncy( 1)=   23;
title_y( 1,:)='Water sound speed (m/s)                 ';
 nplots=0;
   npdpsamp( 1)=        51;
 nparm=   4;
 f_min( 1)=   0.149750E+04;
 f_max( 1)=   0.150250E+04;
par2phy( 1)=  2;
   npdpsamp( 2)=        51;
 f_min( 2)=   0.147750E+04;
 f_max( 2)=   0.148250E+04;
par2phy( 2)=  2;
   npdpsamp( 3)=        51;
 f_min( 3)=   0.500000E+01;
 f_max( 3)=   0.100000E+02;
par2phy( 3)=  9;
   npdpsamp( 4)=        51;
 f_min( 4)=   0.100000E-01;
 f_max( 4)=   0.100000E+03;
par2phy( 4)=  8;
  res=[ 
       1499.500          0.113
       1481.700          0.115
       9300.000          0.000
         78.002          0.000
  ]; 
  bestfit=   1.652139E-04 ;
