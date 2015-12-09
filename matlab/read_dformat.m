function [pres] = read_dformat(p1,p2);
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

[fid,msg] = fopen(file,'r');
if(fid == -1) display(['The file ' file ' could not be opened']);
 end
%keyboard

okay  = 0;   % end-of-file flag

%  -----------------------------------------------------------
%  R E A D I N G   T H E   D A T A   F R O M   T H E   F I L E
%  -----------------------------------------------------------
s = fgetl(fid);
%  s = fgetl(fid); fprintf(1,'%s\n',s);
pres=[]
while (s(1)~=-1) 
  while (s(1)=='!') 
    s = fgetl(fid); fprintf(1,'%s\n',s);
  end
  aux= sscanf(s,' %d %d %d ( %f, %f) ');
  if (aux(5)~=0.)
    pres=[pres aux(4)];
  else
    pres=[pres aux(4)+j* aux(5)];
  end
    s = fgetl(fid);
end
%  for nn=1:Nsen
%        s = fgetl(fid);
%        aux = sscanf(s,' %d (%f,%f) %f');
%        if(nn ~= aux(1) ) nn,aux,keyboard; end;
%keyboard
%        Pvec(nn,ioffset) = aux(2) + j*aux(3);
%     end;
%  end;
  fprintf(2,'whole pressure vector read\n');
%end;






