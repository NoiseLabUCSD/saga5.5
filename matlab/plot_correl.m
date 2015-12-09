%  PLOT_CORREL plots the correlation for the estimated parameters from
%              SAGA. It reads from the file *.cor 
%              There is no arguments to PLOT_CORREL.

if (exist('inputfil')==0)
 filenm=input('SAGA filename ? ','s')
else
 filenm=inputfil
end;
 fil2  =[filenm,'.cor'];
% fil2 ='fort.15'
 fid2  = fopen(fil2 ,'r'); 
 if fid2==-1,
   x=['The file ',fil2,' does not excist']
 end
 clear cov
     cov = fscanf(fid2,'%f',  [nparm,nparm]); 

cov(:,nparm+1)=cov(:,nparm);
cov(nparm+1,:)=cov(nparm,:);
xx=[0:nparm] +0.5;
figure(3)
pcolor(xx,xx,cov)
%colormap('jet')
colormap('gray')
colorbar
xlabel('parameter')
ylabel('parameter')
title(filenm)
shading flat
   set(gca,...
        'ytick',[1:nparm],... 
        'xtick',[1:nparm]);
