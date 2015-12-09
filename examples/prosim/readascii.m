function [data,par]=readtrf(filename);

%--------------------
% USE trftoascii
% TO CONVERT FROM TRF
% TO ASCII FIRST!
%--------------------

disp(' -------------- ');
disp(['reading ' filename]);
disp(' -------------- ');

fid	= fopen([filename]);
if fid<1, error('   '); end;

%--------------
% Read Headers
%--------------


fileid	= fscanf(fid,'%s',1);
disp(fileid)
prognm	= fscanf(fid,'%s',1);
disp(prognm);

nout	= fscanf(fid,'%d',1);
nout
iparm	= fscanf(fid,'%d',nout);
iparm
title	= fscanf(fid,'%c',64);
%title   = fgetl(fid); 
disp(title);

signn	= fscanf(fid,'%s',1);

signn

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

disp(dt)

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


data	= fscanf(fid,'%f',[2*nout,(mx-lx+1)*msuft*nplots*ir]);
fclose(fid);

data	=data';

%s	= data(:,2:2:2*nout).*exp(i*data(:,1:2:2*nout-1));
s	= data(:,1:2:2*nout-1)+1i*data(:,2:2:2*nout);
clear data;
%
data=reshape(s,nplots*ir,mx-lx+1)';
%
par	= [nout];
par	= [par  rd rdlow ir];
par	= [par  r0 rspace nplots];
par	= [par  nx lx mx dt];
par	= [par  msuft];
%par
