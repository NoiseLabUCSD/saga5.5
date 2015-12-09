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
nhori=nplot_c;
nvert=1;
if (nhori>3)
   nvert=2; nhori=(nhori+1)/2;
end;
for icont=1:nplot_c
    subplot(nhori,nvert,icont)
%
%   read sector
%

    sector = fscanf(fid,'%f',  28); % sector contains 28 points
%
%
    nx=npx(icont);
    ny=npy(icont);
    x=eval(xp( icont,:));
    y=eval(yp( icont,:));
    c_val = fscanf(fid,'%f',  [nx,ny]);
%    if (icont==1)   % to make reverse color scale		 
%       c_val=c_val;
%    end
   cmin=-2;
   c_val(c_val<cmin)=ones(size(c_val(c_val<cmin)))*cmin;   
    imagesc(x,y,rot90(c_val(:,:))),colormap(jet),shading('flat'),h=colorbar; 
%   pcolor(x,y,rot90(c_val(:,:))),colormap(jet), shading('interp'),colorbar 
    xlabel(title_x( icont,1:ncx(icont)))
    ylabel(title_y( icont,1:ncy(icont)))
    pos=get(h,'position'); pos(3)=0.015; pos(1)=pos(1)-0.006;
    set(h,'position',pos,'xlim',[0 1])
    set(gca,'Ydir','rev')
end; %each plot

%if(exist('stringofcomments') == 0)
%   stringofcomments0 = ['Results from data ',filenm];
%else
%   stringofcomments0 =  stringofcomments
%end;

%if(exist('no_title') == 0)
% figure(1)
%   axes('position',[0,0,1,1]); axis('off')
%   text(.05,.05, [stringofcomments0],'fontsize',10)
%   text(.05,.02, ['File: ',pwd,'/',filenm,'; Date: ',date],'fontsize',6)
%end;

%end;  %contsaga
orient('tall')





