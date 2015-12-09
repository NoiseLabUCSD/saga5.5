function [pres,par] = read_it(p1,p2);
%[Pvec,freq,zvec,file] = read_pvec(file,Nfreq);
%
% Reading the pressure vector  from a file written in Saga-format d
%
% SYNOPSIS:
%   read file 'pres.in':        [..] = covread(Nfreq);   { default action }
%   read file 'myfile.in'       [..] = covread('myfile.in',Nfreq);
%   
% output parameters:
%   Pvec    =  Ndep X Nfreq   pressure vectors
%   freq    =  1    x Nfreq   frequency
%   zvec    =  1    x Ndep    sensor depths in meters
%   file    =  Filename
%
%  Christoph Mecklenbrauker Oct 95
%  Peter Gerstoft Mar 97

if(nargin==0) % no input arguments: setup defaults
   file = 'pres.in'; 
else
   file = p1; 
 end
igzip=0;
names=dir([file '*']);
%if (exist('names.name')==0)
%  display(['no such file exist ' file]) 
%  keyboard
%end
%if ~isempty(findstr(names.name,[file '.gz']));
if ~isempty(findstr(names.name,['.gz']));
  igzip=1;
  eval(['!gunzip ' file]);
end

%keyboard
[fid,msg] = fopen(file,'r');
if(fid == -1)  error(msg); end

okay  = 0;   % end-of-file flag

%  -----------------------------------------------------------
%  R E A D I N G   T H E   D A T A   F R O M   T H E   F I L E
%  -----------------------------------------------------------
s = fgetl(fid);
%  s = fgetl(fid); fprintf(1,'%s\n',s);
freq  = fscanf(fid,'%f',1);
rmin  = fscanf(fid,'%f',1);
drout = fscanf(fid,'%f',1);
nr    = fscanf(fid,'%d',1);
zrmin = fscanf(fid,'%f',1);
zrinc = fscanf(fid,'%f',1);
nrdep = fscanf(fid,'%d',1);

%keyboard
pres=[];
%scanf(fid,' %d %d %d  %d %d %d %d %d %d    ',par]

data    = fscanf(fid,'%f',[nrdep*2,nr]);
pres       = data(1:2:nrdep*2-1,:)+i*data(2:2:nrdep*2,:);

par=[freq rmin drout nr zrmin zrinc nrdep ];

%fprintf(2,'whole pressure vector read\n');
fclose(fid);

if igzip==1
  eval(['! gzip ' file ' & ']);
end





