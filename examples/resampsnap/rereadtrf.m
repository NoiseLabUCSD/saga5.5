function [alldata,par]=rereadtrf(filename);

% par contains the parameters used in the processing
% par(1) = 1
% par(2) = first receiver height
% par(3) = last receiver height
% par(4) = number of receivers in height
% par(5) = first range
% par(6) = delta range
% par(7) = number of ranges
% par(8) = nx, number of time samples, T=nx*(delta t)
% par(9) = lx, first sample computed by tpem
% par(10) = mx, last sample computed by tpem
% par(11) = delta t
% par(12) = 1
%%--------------------
% USE trftoascii
% TO CONVERT FROM TRF
% TO ASCII FIRST!
%--------------------

disp(' -------------- ');
disp(['reading ' filename]);
%disp(' -------------- ');

fid	= fopen([filename]);
if fid<1, error('   '); end;

%--------------
% Read Headers
%--------------


[fileid,dd]	= fscanf(fid,'%s',1);
%disp(fileid)
prognm	= fscanf(fid,'%s',1);
%disp(prognm);

nout	= fscanf(fid,'%d',1);

iparm	= fscanf(fid,'%d',nout);

title	= fscanf(fid,'%c',64);
%title   = fgetl(fid); 
%disp(title);

signn	= fscanf(fid,'%s',1);


fctrf	= fscanf(fid,'%f',1) ;

sd	= fscanf(fid,'%f',1);
rd	= fscanf(fid,'%f',1);
rdlow	= fscanf(fid,'%f',1);

ir	= fscanf(fid,'%d',1);

% if ir < 0 nog 3? getallen lezen

r0	= fscanf(fid,'%f',1);
rspace	= fscanf(fid,'%f',1);
nplots	= fscanf(fid,'%d',1);

nx	= fscanf(fid,'%d',1);
lx	= fscanf(fid,'%d',1);
mx	= fscanf(fid,'%d',1);
dt	= fscanf(fid,'%f',1);

%disp(dt)

icdr	= fscanf(fid,'%d',1);
omegim	= fscanf(fid,'%f',1);

msuft	= fscanf(fid,'%d',1);
isrow	= fscanf(fid,'%d',1);
inttyp	= fscanf(fid,'%d',1);

nopp	= fscanf(fid,'%d',2);
nopp	= fscanf(fid,'%f',5);

%-------------------
% read data
%-------------------

%
par	= [nout];
par	= [par  rd rdlow ir];
par	= [par  r0 rspace nplots];
par	= [par  nx lx mx dt];
par	= [par  msuft];
%par
filenamebin='fort.34'
fidb     = fopen([filenamebin]);
if fidb<1, error(['file ',filenamebin,' could not be opened   ']); end;
[s,w] = unix('wc fort.81')
ntrf=str2num(w(1:8))
datasize=ir*nplots*(mx-lx+1)*2; %for complex
nfreq=(mx-lx+1);
alldata=zeros(ir,nplots,ntrf);

%keyboard  
backspacestring = 8 * ones(1,40);

for ii=1:ntrf
  fprintf(1,'Reading trf number = %d',ii)
  dum   = fread(fidb,[ 1],'float');
  stemp   = fread(fidb,[ ir*2,nplots],'float');
  dum   = fread(fidb,[ 1],'float');
  stemp_c= stemp(1:2:end,:)+i*stemp(2:2:end,:);
% stemp_tl=20*log10(abs(stemp_c));
  alldata(:,:,ii)=stemp_c;
  fprintf(1,'%s',backspacestring);
end
  
%if (rem(ii,100)==1); fprintf(1,' %d',ii); end;
%%%%%  
%[fileid,dd]	= fscanf(fid,'%s',1);
%[data, dd]	= fscanf(fid,'%f',[datasize]);

%if (dd==0)
%  disp(['read whole trf file, with ',int2str(ii-1),' environments'])
%  return; end

%data=reshape(data,[ir,nplots,nfreq]);
%data_r=reshape(data(1:2:end),[ir,nplots,1]);
%data_i=reshape(data(2:2:end),[ir,nplots,1]);
%data_c=data_r+i*data_i;
%data_tl=20*log10(abs(data_c))
%alldata(:,:,ii)=stemp_c;
%    backspacestring = 8 * ones(1,40);
%    fprintf(1,'%s',backspacestring);
 
%end









