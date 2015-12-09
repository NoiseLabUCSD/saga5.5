function [ranges,pres,comment] = read_Dformat(p1);
%[ranges,pres,comment] = read_Dformat(p1);
%
% Reading the pressure vector  from a file written in Saga-format d
%
% SYNOPSIS:
%   read file 'pres.in':        [..] = covread(Nfreq);   { default action }
%   read file 'myfile.in'       [..] = covread('myfile.in');
%   
% output parameters:
%   Pvec    =  Nrange   pressure vectors
%   Pvec    =  Nrange   pressure vectors
%
%  Peter Gerstoft Sep 00
if(nargin==0) % no input arguments: setup defaults
   file = 'pres.in'; 
else
   file = p1; 
end

[fid,msg] = fopen(file,'r');
if(fid == -1)  display(['The file ' file ' could not be opened']); end

okay  = 0;   % end-of-file flag

%  -----------------------------------------------------------
%  R E A D I N G   T H E   D A T A   F R O M   T H E   F I L E
%  -----------------------------------------------------------
s = fgetl(fid);
%  s = fgetl(fid); fprintf(1,'%s\n',s);
pres=[]; ranges=[]; comment=[]; 
while (s(1)~=-1) 
  while (s(1)=='!') 
    comment=[comment s]
    s = fgetl(fid); fprintf(1,'%s\n',s);
  end
  if (s(1)~=-1) 
    aux= sscanf(s,'%f  %f ');
    ranges=[ranges aux(1)];
    pres=[pres aux(2)];
    s = fgetl(fid);
  end
end

  fprintf(2,'whole pressure vector read\n');
%end;






