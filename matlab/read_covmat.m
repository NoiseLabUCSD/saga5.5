function [Cov,freq,zvec,covfile] = read_covmat(p1,p2);
%[Cov,freq,zvec,covfile] = read_covmat(covfile,Nfreq);
%
% Reading the covariance matrix in from a file written in Saga-format
%
% SYNOPSIS:
%   read file 'cov.in':        [..] = read_covmat;   { default action }
%   read file 'cov_dpss.in'    [..] = read_covmat('./ElbaBB_data/cov_dpss.in',3);
%   
% output parameters:
%   Cov     =  Ndep x (Ndep*Nfreq)  concatenated estimated covariance matrices
%   freq    =  1  x Nfreq	    frequency
%   zvec    =  1  x Ndep	    sensor depths in meters
%   pwrl    =  1  x Nfreq	    power levels (not used)
%   covfile =  Filename
%
%  Christoph Mecklenbrauker Oct 95
%  Peter Gerstoft Mar 97
if(nargin==0) % no input arguments: setup defaults
   covfile = './ElbaBB_data/cov.in'; Nfreq = 6;
   cleanup = 1; method = 2;
elseif(nargin==2)
   covfile = p1; Nfreq = p2; cleanup = 1; method = 2;
else
   covfile = p1; Nfreq = p2; cleanup = p3; method = p4;
end

[fid,msg] = fopen(covfile,'r');
if(fid == -1)  display(['The file ' file ' could not be opened']); end

%Nsen  = 48;  % Number of sensors in array
 aux(1:2)=[0 0]
okay  = 0;   % end-of-file flag


%  -----------------------------------------------------------
%  R E A D I N G   T H E   D A T A   F R O M   T H E   F I L E
%  -----------------------------------------------------------
for ifreq=1:Nfreq
  s = fgetl(fid); fprintf(1,'Title1: %s\n',s);
  if(strcmp(s(1:32),'! File format: Covariance matrix')==1)
   s = fgetl(fid);
   while (s(1)=='!')
     s = fgetl(fid); fprintf(1,'%s\n',s);
   end  
   fprintf(1,'Title2: %s\n',s);
 else
%%%%%    fprintf('
  end
%  if(strcmp(s(1:8),'26-16b2f')==0) error('wrong file'); end;
  s = fgetl(fid);
  aux = sscanf(s,'%f %f');
  freq(ifreq) = aux(1);  fprintf(1,'Frequency: %f\n', freq(ifreq));
%  pwrl(ifreq) = aux(2);
   nvec=1
%   if (aux(2)>0 ) nvec=aux(2), end;
  s = fgetl(fid);
  aux = sscanf(s,'%f %f');
  aux
  Nsen = aux(1) ;  fprintf(1,' No of sensors: %f\n', Nsen);
     if (size(aux,2)>1 ) Mfreq=aux(2); 
		  fprintf(1,' No of frequencies in file: %f\n',Mfreq); 
     end;
%  n = aux(1); if(n ~= Nsen) s,n,error('n ~= 48'); end;
  zvec = fscanf(fid,' %f',Nsen);
  s = fgetl(fid); % read the EOL character properly
  offset = (ifreq-1)*Nsen;
%
%  Now the covariance matrix is finally read 
%
  for nn=1:Nsen*nvec
     fprintf(2,'%3d\b\b\b',nn);
     for mm=1:Nsen*nvec
        s = fgetl(fid);
        aux = sscanf(s,' %d %d (%f,%f)');
        if(nn ~= aux(1) | mm ~= aux(2)) mm,nn,aux,keyboard; end;
%        Cov(nn,mm+offset) = aux(3) + j*aux(4);
        Cov(nn,mm,ifreq) = aux(3) + j*aux(4);
     end;
  end;
  fprintf(2,'whole covariance matrix read\n');
end;






