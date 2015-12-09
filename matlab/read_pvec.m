function [Pvec,freq,zvec,file] = read_pvec(p1,p2);
%[Pvec,freq,zvec,file] = read_pvec(file,Nfreq);
%
% Reading the pressure vector  from a file written in Saga-format
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
if(nargin==1) % no input arguments: setup defaults
   file = 'pres.in'; Nfreq = p1;
else
   file = p1; Nfreq = p2;
end

[fid,msg] = fopen(file,'r');
if(fid == -1)  display(['The file ' file ' could not be opened']); end

okay  = 0;   % end-of-file flag

%  -----------------------------------------------------------
%  R E A D I N G   T H E   D A T A   F R O M   T H E   F I L E
%  -----------------------------------------------------------

s = fgetl(fid);
while (s(1)=='!')
  s = fgetl(fid); fprintf(1,'%s\n',s);
end  
for ifreq=1:Nfreq
  s = fgetl(fid);

  aux = sscanf(s,' %f ');
  freq(ifreq) = aux(1);
  s = fgetl(fid);
  Nsen = sscanf(s,'%d'); 
  fprintf(1,'Frequency %d Sensors %d \n',freq(ifreq),Nsen);
  zvec = fscanf(fid,' %f',Nsen);
  s = fgetl(fid); % read the EOL character properly
  ioffset = ifreq;
  for nn=1:Nsen
    s = fgetl(fid);
    aux = sscanf(s,' %d (%f,%f) %f');
    if(nn ~= aux(1) ); nn,aux,keyboard; end;
    %keyboard
    Pvec(nn,ioffset) = aux(2) + j*aux(3);
  end;
end;
fprintf(2,'whole pressure vector read\n');







