% CONTSAGA plots the contour output from the SAGA postprocessor.
% 
% contsaga will inquire about the file from which to plot the data.
%          alternatively, the filename can be specified in the
%          variable 'inputfil'.
%          By giving the variable "no_title" a value, title will be
%          plotted at the bottom. 

%function plotsaga_1(inputfil,stringofcomments)
% filenm='ys3';
if (exist('inputfil')==0)
   filenm=input('SAGA filename ? ','s')
else
   filenm=inputfil
end;
fil  =[filenm,'.bdr'];
fid  = fopen(fil ,'r');
if fid==-1,
     x=['The file ',fil,' does not excist']
end
%
 clear c_val ; 
%
eval(filenm) 
iforward= iopt(30);

for icont=1:nplot_c
subplot(nplot_c,1,icont)
%
% read sector
%

     sector = fscanf(fid,'%f',  28); % sector contains 28 points
%
%
nx=npx(icont)
ny=npy(icont)
x=eval(xp( icont,:));
y=eval(yp( icont,:));
     c_val = fscanf(fid,'%f',  [nx,ny]);

     pcolor(x,y,rot90(c_val(:,:))),colormap(gray), shading('interp'),colorbar 
%     pcolor(x,y,rot90(c_val(:,:))),colormap(jet), shading('interp'),colorbar 
xlabel(title_x( icont,:))
ylabel(title_y( icont,:))
end; %each plot

if(exist('stringofcomments') == 0)
 stringofcomments = ['Results from data ',filenm];
end;
if(exist('no_title') == 0)
 figure(2)
 axes('position',[0,0,1,1]); axis('off')
 text(.05,.05, [stringofcomments],'fontsize',10)
 text(.05,.02, ['File: ',pwd,'/',filenm,'; Date: ',date],'fontsize',6)
 figure(1)
 axes('position',[0,0,1,1]); axis('off')
 text(.05,.05, [stringofcomments],'fontsize',10)
 text(.05,.02, ['File: ',pwd,'/',filenm,'; Date: ',date],'fontsize',6)
end;

end;






