% PLOTSAGA plots the output from the SAGA postprocessor.
% 
% Plotsaga will inquire about the file from which to plot the data.
%          alternatively, the filename can be specified in the
%          variable 'inputfil'.
%          By giving the variable "no_title" a value, title will be
%          plotted at the bottom. 

%function plotsaga_1(inputfil,stringofcomments)
% inputfil='ys3';
clear all;
if (exist('inputfil')==0)
   filenm=input('SAGA filename ? ','s')
else
   filenm=inputfil;
end;
fil  =[filenm,'.plt'];
fid  = fopen(fil ,'r');
if fid==-1,
     x=['The file ',fil,' does not exsist']
end
%
 clear pdp a plaxis pdpref pdpxxx pdpref1 pdpxxx1 ; 
%
unitplot=0;
eval(filenm)
iforward= iopt(30);
loneplot = iopt(14); 
if (unitplot==1),
   loneplot =0;
end
%
% read fitness plots
%
if iopt(10)>0
%   if iopt(30)==5
%      a = fscanf(fid,'%f',  [nx,nplots*2*2]); %2*2: 2lineplots xand y val.
%   else
      a = fscanf(fid,'%f',  [nphone,nplots*2*2]); %2*2: 2lineplots xand y val.
%   end 
end
%
%
if loneplot==0,                   % condensed plot
    dx       = 1/(npdpsamp-1); 
    xtexts   = 'range of parameter';
    pdp      = fscanf(fid,'%f',[npdpsamp(1),nparm]); 
    pdp      = pdp';
    if iopt(19)==1,        % option A
        pdpref   = fscanf(fid,'%f',[2*nparm]);   % These are the ref values
        pdpxxx   = fscanf(fid,'%f',[2*nparm]);   % for option A
    end;
else                              % one plot for each parameter
    pdp=zeros(nparm,max(npdpsamp));
    for j=1:nparm
       pdp(j,1:npdpsamp(j))    = fscanf(fid,'%f',[1,npdpsamp(j)]);
       if iopt(19)==1,        % option A
          pdpref1(:,j)= fscanf(fid,'%f',[2]);
          pdpxxx1(:,j)= fscanf(fid,'%f',[2]);
       end;
    end 
  end;
%%%%%%%%%%
%%%%%
%%%%%  DO not plot anything if noplot excist
%%%%%
  
fclose(fid);
if (exist('noplot')==0)
%
% plot pdp
%
 figure(1)
 clf
 if (loneplot==0)                    % condensed plot
   if (unitplot==0) 
       hwiggle2(pdp,[dx	1], [0 1],2,'b')
       axis ([ 0 1 0 nparm])
       xlabel(xtexts)
   else
     dx=(f_max(1)-f_min(1))/(npdpsamp(1)-1);
     test2=f_min(1):dx:f_max(1);
       hwiggle2(pdp,[dx	1], [f_min(1) 1],2,'b')
     axis ([ f_min(1) f_max(1)  0 nparm+1])
     xlabel( xtitles(iforward,par2phy(1)),'Fontsize',12);
   end
   hold;
   if iopt(19)==1,        % option A
       plot(pdpref,pdpxxx,'r');
   end;
   hold;
else                               % one plot for each parameter
 
  left=.18;
  bottom0=1.0;
  width=.6;
  height=.15;
  fact=2.6;
  xhelp=fact*height;

  if (xhelp*nparm >.92)
      height=.92/fact/nparm
  end

  for j=1:nparm
     bottom=bottom0-(j)*(height*fact);
     subplot('position',[left bottom width height]);
%%%%%     subplot(2*nparm,1,2*j)
     dx=(f_max(j)-f_min(j))/(npdpsamp(j)-1);
     if (0==1)
       dx=dx*2;
       testppd=pdp(j,1:2:npdpsamp(j))+pdp(j,2:2:npdpsamp(j));
     else
       testppd=pdp(j,1:npdpsamp(j));       
     end
     test2=f_min(j):dx:f_max(j);
     wplot(test2,testppd,'b');
     axis ([ f_min(j) f_max(j)  0 1])
     xlabel( xtitles(iforward,par2phy(j)),'Fontsize',10);
     set(gca,'Ytick',2);
     hold on;
     if iopt(19)==1,        % option A
         plot(pdpref1(:,j),pdpxxx1(:,j),'r');
     end;
     hold off;
     set(gca,'Fontsize',11);
   end
 end
 if (nparm>6 & unitplot==0)
   orient tall
 end
%
%
%
 if nplots>0
 mhoriz=3;
 mvert=ceil(nplots/mhoriz);
 if (nplots < mhoriz) mhoriz=nplots; end;
 nstep_plots=1;
 if (nplots>20)
    nstep_plots=round(nplots/20);
    disp(['Only every ',int2str(nstep_plots),' curve is plotted'])
 end
 mvert=ceil(nplots/mhoriz/nstep_plots);
 
 figure(2)
 clf
 jplot=0;
% keyboard
 for j=1:nstep_plots:nplots
   jplot=jplot+1;
   subplot(mvert,mhoriz,jplot)
%%%%%%%%%%% To plot more squeezed plots
%   apos=get(gca,'position');
%  set(gca,'position',[apos(1)*.95 apos(2)*.95 apos(3) apos(4)]);
%%%%%%%%%%%
   j1=4*(j-1);

   if (iforward==1 | iforward==2 | iforward==3| (iforward==7 & ...
	 iopt(5)==4) | (iforward==8 & iopt(5)==4)| iforward==9 ),
%    z=1:32;
%   plot(a(:,j1+1),z,'y-',a(:,j1+3),z,'m-.')
%   axis([-0.0 0.1 1 32]);
     plot(a(:,j1+1),a(:,j1+2),'g-',a(:,j1+3),a(:,j1+4),'b--')
     if (plaxis(4,j)< plaxis(3,j));
        dum= plaxis(4,j);
        plaxis(4,j)= plaxis(3,j);
        plaxis(3,j)= dum;
        set(gca,'Ydir','reverse');
      end
%      keyboard
  elseif ( (iforward==7 | iforward==8) &   iopt(5)==5),  % horizontal array
     plot(a(:,j1+1),a(:,j1+2),'g-',a(:,j1+3),a(:,j1+4),'b--')
  elseif (iforward==4| iforward==7 | iforward==8 ),
     plot(a(:,j1+1),a(:,j1+2),'g-',a(:,j1+3),a(:,j1+4),'b--')
  elseif iforward==6,    %TPEM
     if iopt(3)==12,
      plot(a(:,j1+2),a(:,j1+3),'g-',a(:,j1+4),a(:,j1+3),'b--')
      xtemp=plaxis(:,j);
      plaxis(1,j)=xtemp(3);
      plaxis(2,j)=xtemp(4);
      plaxis(3,j)=xtemp(1);
      plaxis(4,j)=xtemp(2);       
     else
      plot(a(:,j1+1),a(:,j1+2),'g-',a(:,j1+3),a(:,j1+4),'b--')
        dum= plaxis(4,j);
        plaxis(4,j)= plaxis(3,j);
        plaxis(3,j)= dum;
%        plaxis(4,j)= plaxis(3,j);
%        plaxis(3,j)= 0;                %plaxis(4,j)
     end
   elseif iforward==5,
     plot(a(:,j1+1),a(:,j1+2),'g-',a(:,j1+3),a(:,j1+4),'b--')
   end
  xrev=0;  yrev=0;
  if  (plaxis(4,j)< plaxis(3,j))
      xtemp=plaxis(:,j);
      plaxis(3,j)=xtemp(4);
      plaxis(4,j)=xtemp(3);
      yrev=1
  end
  if  (plaxis(2,j)< plaxis(1,j))
      xtemp=plaxis(:,j);
      plaxis(1,j)=xtemp(2);
      plaxis(2,j)=xtemp(1);
      xrev=1
  end
    axis(plaxis(:,j));
  if (xrev==1) set(gca,'xdir','rev'); end
  if (yrev==1) set(gca,'ydir','rev'); end
   xtext=['curve ',num2str(j)];
    set(gca,'Fontsize',11);
%   text(0.1,32,xtext,'VerticalAlignment','top','HorizontalAlignment','r')
 end
 if (mvert>3) 
   orient tall;
 end
end
%%%
%  Put comments and date on each figure:
%%%
if(exist('stringofcomments') == 0);
 stringofcomments0 = ['Results from data ',filenm];
else
 stringofcomments0 =  stringofcomments;
end;

%nhq 20140902
% if(exist('no_title') == 0)
%  figure(2)
%  axes('position',[0,0,1,1]); axis('off');
%  text(.05,.05, [stringofcomments0],'fontsize',10);
%  text(.05,.02, ['File: ',pwd,'/',filenm,'; Date: ',date],'fontsize',6);
%  figure(1)
%  axes('position',[0,0,1,1]); axis('off');
%  text(.05,.05, [stringofcomments0],'fontsize',10);
%  text(.05,.02, ['File: ',pwd,'/',filenm,'; Date: ',date],'fontsize',6);
% end;
% 
end
