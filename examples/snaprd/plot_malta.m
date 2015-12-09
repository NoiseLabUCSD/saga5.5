% PLOTSAGA plots the output from the SAGA postprocessor.
% 
% Plotsaga will inquire about the file from which to plot the data.
%          alternatively, the filename can be specified in the
%          variable 'inputfil'.
 

%function plotsaga_1(inputfil,stringofcomments)
 filenm='tl_malta_rd';
%if (exist('inputfil')==0)
%   filenm=input('SAGA filename ? ','s')
%else
%   filenm=inputfil
%end;
fil  =[filenm,'.plt'];
fid  = fopen(fil ,'r');
if fid==-1,
     x=['The file ',fil,' does not excist']
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
if iopt(10)==1
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
%
%
%
 if nplots>0
 mhoriz=4;
 mvert=ceil(nplots/mhoriz);
 if (nplots < mhoriz) mhoriz=nplots; end;
 nstep_plots=1
 if (nplots>20)
    nstep_plots=round(nplots/20);
    disp(['Only every ',int2str(nstep_plots),' curve is plotted'])
 end
 mvert=ceil(nplots/mhoriz/nstep_plots)

 figure(2)
 clf
 jplot=0;
 for j=1:nstep_plots:nplots
   jplot=jplot+1;
   subplot(mvert,mhoriz,jplot)
  apos=get(gca,'position');
  set(gca,'position',[apos(1)*.95 apos(2)*.95 apos(3) apos(4)]);
   j1=4*(j-1);
   if (iforward==1 | iforward==2 | iforward==3| (iforward==7 & iopt(5)==4) ),
%    z=1:32;
%   plot(a(:,j1+1),z,'y-',a(:,j1+3),z,'m-.')
%   axis([-0.0 0.1 1 32]);
     plot(a(:,j1+1),a(:,j1+2),'g-',a(:,j1+3),a(:,j1+4),'b-')
     if (plaxis(4,j)< plaxis(3,j)),
        dum= plaxis(4,j);
        plaxis(4,j)= plaxis(3,j);
        plaxis(3,j)= dum;
        set(gca,'Ydir','reverse');
     end
  elseif (iforward==4| iforward==7 ),
     plot(a(:,j1+1),a(:,j1+2),'g-',a(:,j1+3),a(:,j1+4),'m--')
  elseif iforward==6,
     plot(a(:,j1+1),a(:,j1+2),'g-',a(:,j1+3),a(:,j1+4),'m--')
     plaxis(4,j)= plaxis(3,j);
     plaxis(3,j)= 0;                %plaxis(4,j)
   elseif iforward==5,
     plot(a(:,j1+1),a(:,j1+2),'g-',a(:,j1+3),a(:,j1+4),'m--')
   end
   axis(plaxis(:,j));
   axis([0 17 45 81])
   if (j<nplots-mhoriz+1)
        set(gca, 'xtick',[]);
   end
   if (~(rem(j,mhoriz)==1))
        set(gca, 'ytick',[]);
   end
   xtext=['curve ',num2str(j)];
     set(gca,'Fontsize',9);
%   text(0.1,32,xtext,'VerticalAlignment','top','HorizontalAlignment','r')
 end
 if (mvert>3) 
 orient tall
 end
end
%%%
%  Put comments and date on each figure:
%%%
if(exist('stringofcomments') == 0)
 stringofcomments = ['Results from data ',filenm];
end;
 figure(2)
 axes('position',[0,0,1,1]); axis('off')
 text(.45,.05, 'Range (km)','fontsize',13)
 text(.08,.35, 'Transmission loss (dB)','Rotation',90,'fontsize',13)
% text(.05,.02, ['File: ',pwd,'/',filenm,'; Date: ',date],'fontsize',6)
















