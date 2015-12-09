function  [M,z]= trilin(baseheight,thick,Mdef,c1,delta)
%function trilin(baseheight,thick,Mdef,c1)
% refractivity profile 
% Peter Gerstoft, Aug 00
%
%these 8 parameters specify the profile
c0=0.13;         % slope in evaporation and mixed layer adiabatic
%c1=0         % slope in evaporation and mixed layer adiabatic
c2= 0.118;       % slope above inversion layer
%delta=50        % evaporation duct
z0=0.00015;      % roughness
%baseheight=100  %
%%%%%thick=50
offset=330;
%Mdef= 10
%these are specified by the PE
zmax=500;      % maximum heigth of PE
dz=0.7;  % pe discretization



%evaporation duct
ahelp=1/(1-c1/c0);
if (ahelp>0 & ahelp<2) 
  zd=ahelp*delta;
else
  zd=2*delta;
end

zd=min(zd,baseheight);
if (zd>dz) 
  z=0:dz:(zd+dz/2); z(1)=0.01;
  M=c0*(z-zd-delta*log(z/zd));
else
  zd=0;
  z(1)=0;
  M(1)=0;
end;
iz=floor((zd+dz/2)/dz+1);
Mr = M(iz); zr = z(iz);

%keyboard
%mixed layer
baseheight_trunc=max(baseheight,0);
if (baseheight>0) 
  zloc =(zr-zd+dz):dz:(baseheight-zd+dz/2);
  z=[z zr+zloc];
  M=[M c1*zloc+Mr];
  iz=floor((baseheight+dz/2)/dz+1);
end  
Mr = M(iz); zr = z(iz);

% inversion layer
thick_trunc=max(thick,0);
if (thick>0) 
  zh = baseheight_trunc+thick+dz/2;  %top in global coordinates
  zloc =(zr-baseheight_trunc+dz):dz:(thick+dz/2); 
  z =[z zr+zloc];
  dM=-Mdef/thick;
  M=[M  Mr+dM*zloc];
  izh=floor(((thick+dz/2)-(zr-baseheight_trunc))/dz);
else
  izh=0  
end

iz2=iz+izh;
Mr=M(iz2); zr = z(iz2);

% top layer
zup=zmax-baseheight_trunc-thick_trunc -zd;
zloc =dz:dz:zup;
z =[z zr+zloc];
M=[M  Mr+c2*zloc];

% Galinas offset
%offset2=offset-c1*(delta-delta*log(delta/z0));

% Is this the definition of offset
%offset2=offset+c1*(delta*log(delta/z0));

M=M+offset;
if nargout == 0
plot(M,z)
ylabel('Height (m)')
xlabel('Modified refractivity (M-units)')
end






