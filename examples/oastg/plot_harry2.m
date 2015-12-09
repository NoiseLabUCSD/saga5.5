% CONTSAGA plots the contour output from the SAGA postprocessor.
% 
% contsaga will inquire about the file from which to plot the data.
%          alternatively, the filename can be specified in the
%          variable 'inputfil'.
%          By giving the variable "no_title" a value, title will be
%          plotted at the bottom. 

%function plotsaga_1(inputfil,stringofcomments)
%inputfil='varray_10db';
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
nx=npx(icont);
ny=npy(icont);
x=eval(xp( icont,:));
y=eval(yp( icont,:));
     c_val = fscanf(fid,'%f',  [nx,ny]);
cm=30;
if (icont==2) c_val(c_val>cm)=ones(size(c_val(c_val>cm)))*cm;   end;
if (icont==3)  c_val(c_val>cm)=ones(size(c_val(c_val>cm)))*cm;   end;
     pcolor(x,y,rot90(c_val(:,:))),colormap(jet), shading('flat')
%     pcolor(x,y,rot90(c_val(:,:))),colormap(jet), shading('interp'),colorbar 
if (icont==1)  
  c_val=-c_val;
  cmax=max(max(c_val));
  c_val= c_val-cmax;
  cmin=-2.5;
   c_val(c_val<cmin)=ones(size(c_val(c_val<cmin)))*cmin;   
  pcolor(x,y,rot90(c_val(:,:))),colormap(jet), shading('flat'),h=colorbar;
end;
if (icont==2 )   
       pcolor(x,y,rot90(c_val(:,:))),colormap(jet), ...
               shading('flat'),caxis([0 cm]);h=colorbar('vert');  end;
if ( icont==3)  
       pcolor(x,y,rot90(c_val(:,:))),colormap(jet), ...
               shading('flat'),caxis([0 cm]);h=colorbar('vert');  end;
 set(gca,'Ydir','rev')
xlabel(title_x( icont,1:ncx(icont)))
ylabel(title_y( icont,1:ncy(icont)))
pos=get(h,'position'); pos(3)=0.025; pos(1)=pos(1)-0.023;
if (icont==2) pos(1)=pos(1)+1; end;
if (icont==3) pos(4)=0.5194; end;
 set(h,'position',pos,'xlim',[0.3 1])
set(gca,'Ydir','rev')
%axis([ 500 2000 0 50])
%keyboard
crb(icont) = mean(mean(c_val)); 
disp('average'), mean(mean(c_val)),
end; %each plot

%CRB=[CRB; crb]

if(exist('stringofcomments') == 0)
 stringofcomments0 = ['Results from data ',filenm];
else
 stringofcomments0 =  stringofcomments
end;
if(exist('no_title') == 0)
 axes('position',[0,0,1,1]); axis('off')
 text(.05,.05, [stringofcomments0],'fontsize',10)
 text(.05,.02, ['File: ',pwd,'/',filenm,'; Date: ',date],'fontsize',6)
end;

orient('tall')

%end;






