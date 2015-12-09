% PLOTSAGA plots the output from the SAGA postprocessor.
% 
% Plotsaga will inquire about the file from which to plot the data.
%          alternatively, the filename can be specified in the
%          variable 'inputfil'.
%          By giving the variable "no_title" a value, title will be
%          plotted at the bottom. 

%function plotsaga_1(inputfil,stringofcomments)
% inputfil='ys3';
if (exist('inputfil')==0)
   filenm = input('SAGA filename ? ','s')
else
   filenm = inputfil;
end;
fil = [filenm,'.plt'];
fid = fopen(fil ,'r');
if fid == -1,
   x = ['The file ',fil,' does not exist']
end 
% unwrap phase using the phone having the maximum energy as reference  
% .cfh.
phaseunwrap = true;
%
 clear pdp a plaxis pdpref pdpxxx pdpopt1 pdpstd1 pdpref1 pdpxxx1 res resloc; 
%
unitplot = 0;
eval(filenm)
temp = res;
if iopt(19)==1,        % option A
   fit=[]; res =[];pointer = [];fitloc=[];resloc=[];
   results;res = [];pointer=[]; fit=[];
   resloc = resloc([end-nparm+1:end]);
end
res = temp;
clear temp
%%%
% set the lables if comparison of data and best model is reguired
%%% iopt(10) > 10
% read fitness plots
if iopt(10) > 0
   a = fscanf(fid,'%f',[nphone, nplots*2*2]); %2*2: 2lineplots xand y val.
   % .cfh.
   if ~exist('freq') 
      freq = input('Signal frequencies: ');
   end    
   nfreq = length(freq);
   for jfre = 1:nfreq
      title_str{jfre} = ['F = ' num2str(freq(jfre)) ' Hz'];
   end
% eliminate the phase difference between data and prediction
% .cfh.
   if phaseunwrap == true;
      iopt_10 = nplots/nfreq; 
      isubplot = [2:iopt_10:nplots];
      tol = 181;  % degree
      for ii = 1:length(isubplot);
         icol = (isubplot(ii)-1)*4+1;
         [ima, iel] = max(a(:,icol-4)); % the ref element is the one that has
                                        % the maximun magnitude
         a(:,icol+2) = a(:,icol+2) - a(iel,icol+2) ; 
         a(:,icol)   = a(:,icol)   - a(iel,icol); 
         resi = a(:,icol)-a(:,icol+2);  %        
         if ~isempty(find(resi > tol))
            a(find(resi > tol),icol+2) = a(find(resi > tol),icol+2)+360;
         end
         if ~isempty(find(resi < -180.1))
            a(find(resi < -tol),icol+2) = a(find(resi < -tol),icol+2)-360;
         end
      end
   end
   jtemp = 0;
   xlabel_str{1} = 'Magnitude';
   xlabel_str{2} = 'Phase (deg)';
   xlabel_str{3} = 'Bartlett Power';
end
iforward = iopt(30);
loneplot = iopt(14); 
if (unitplot==1),
   loneplot = 0;
end
%
if loneplot==0,                   % condensed plot
   dx       = 1/(npdpsamp-1); 
   xtexts   = 'range of parameter';
   pdp      = fscanf(fid,'%f',[npdpsamp(1),nparm]); 
   pdp      = pdp';
   if iopt(19) == 1,        % option A
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
end
% scaling the unit of the SR parameter
% .cfh.
isr = find(par2phy == 9);
temp = repmat(1,nparm,1);
temp(isr) = 0.001;

if iopt(19)==1,        % option A
   pdpstd1 = [res(:,2).*temp - res(:,3).*(f_max-f_min).', res(:,2).*temp  + res(:,3).* ...
              (f_max-f_min).'];
   pdpopt1 = repmat(res(:,1),1,2).';
end
%
% plot pdp
%
figure(1)
clf
if (loneplot == 0)                    % condensed plot
   if (unitplot == 0) 
      hwiggle2(pdp,[dx	1], [0 1],2,'b')
      axis ([ 0 1 0 nparm])
      xlabel(xtexts)
   else
      dx=(f_max(1)-f_min(1))/(npdpsamp(1)-1);
      test2=f_min(1):dx:f_max(1);
      hwiggle2(pdp,[dx	1], [f_min(1) 1],2,'b')
      axis ([ f_min(1) f_max(1)  0 nparm+1])
      xlabel(deblank(xtitles(iforward,par2phy(1))),'Fontsize',10);
   end
   hold;
   if iopt(19)==1,        % option A
      plot(pdpref,pdpxxx,'r');
   end;
   hold;
else                               % one plot for each parameter
   left=.1;
   bottom0=1.0;
   width=.6;
   height=.06;
   fact=2.2;
   xhelp=fact*height;
   if (xhelp*nparm >.92)
      height = .92/fact/nparm;
   end
   
   for j=1:nparm
      bottom = bottom0-(j)*(height*fact);
      subplot('position',[left bottom width height]);
      testppd = pdp(j,1:npdpsamp(j));
%%%%%     subplot(2*nparm,1,2*j)
      dx = (f_max(j)-f_min(j))/(npdpsamp(j)-1);
      test2 = f_min(j):dx:f_max(j);
      wplot(test2,testppd,'b');
      axis ([ f_min(j) f_max(j)  0 1])
      xlabel(deblank(xtitles(iforward,par2phy(j))),'Fontsize',10);
      set(gca,'Ytick',2);
      hold on;

     if iopt(19) == 1,     % option A -- iopt(19) = 1
                           % pdpopt1: SAGA_bestfit model (blue solid)
                           % pdpstd1: one std from the mean (green solid)
                           % res(j,2): SAGA_mean (green solid)
                           % resloc(j): bestfit of Powell method (red solid)
                           % .cfh.
        plot(pdpopt1(:,j).*temp(j),pdpxxx1(:,j),'b','linewidth',1.5);
        plot(pdpstd1(j,:),[mean(pdpxxx1(:,j)) mean(pdpxxx1(:,j))],'-o',...
             'color',[0 0.65 0],'linewidth',1.5,'markersize',4);
        plot(res(j,2).*temp(j),mean(pdpxxx1(:,j)),'+',...
             'color',[0 0.65 0],'linewidth',1.5,'markersize',8)
        plot([resloc(j) resloc(j)].*temp(j),pdpxxx1(:,j),'--r',...
             'linewidth',1.5);
        %        plot(pdpref1(:,j),pdpxxx1(:,j),'--r','linewidth',1.5);
     end;
     hold off;
     set(gca,'Fontsize',10);
   end
end
if (nparm>6 & unitplot==0)
   orient tall
end
%
%
%
if nplots > 0
   mhoriz = 3;
   mvert = ceil(nplots/mhoriz);
   if (nplots < mhoriz) mhoriz=nplots; end;
   nstep_plots = 1;
   if (nplots>20)
      nstep_plots = round(nplots/20);
      disp(['Only every ',int2str(nstep_plots),' curve is plotted'])
   end
   mvert = ceil(nplots/mhoriz/nstep_plots);
   figure(2)
   clf
   jplot = 0;
% keyboard
   for j = 1:nstep_plots:nplots
      jplot = jplot+1;
      subplot(mvert,mhoriz,jplot)
%%%%%%%%%%% To plot more squeezed plots
%   apos = get(gca,'position');
%   set(gca,'position',[apos(1)*.95 apos(2)*.95 apos(3) apos(4)]);
%%%%%%%%%%%
      j1=4*(j-1);

      if (iforward==1 | iforward==2 | iforward==3 | ...
          (iforward==7&iopt(5)==4) | (iforward==8&iopt(5)==4) | ...
          iforward==9 | iforward==11),
         plot(a(:,j1+1),a(:,j1+2),'b-',a(:,j1+3),a(:,j1+4),'r--')
         if (plaxis(4,j) < plaxis(3,j));
            dum = plaxis(4,j);
            plaxis(4,j) = plaxis(3,j);
            plaxis(3,j) = dum;
            %  set(gca,'Ydir','reverse');
         end
         %      keyboard
      elseif ((iforward==7 | iforward==8) & iopt(5)==5),  % horizontal array
         plot(a(:,j1+1),a(:,j1+2),'b-',a(:,j1+3),a(:,j1+4),'r--')
      elseif (iforward==4| iforward==7 | iforward==8 ),
         plot(a(:,j1+1),a(:,j1+2),'b-',a(:,j1+3),a(:,j1+4),'r--')
      elseif iforward==6,    %TPEM
         if iopt(3)==12,
            plot(a(:,j1+2),a(:,j1+3),'b-',a(:,j1+4),a(:,j1+3),'r--')
            xtemp = plaxis(:,j);
            plaxis(1,j) = xtemp(3);
            plaxis(2,j) = xtemp(4);
            plaxis(3,j) = xtemp(1);
            plaxis(4,j) = xtemp(2);       
         else
            plot(a(:,j1+1),a(:,j1+2),'b-',a(:,j1+3),a(:,j1+4),'r--')
            dum = plaxis(4,j);
            plaxis(4,j) = plaxis(3,j);
            plaxis(3,j) = dum;
%        plaxis(4,j)= plaxis(3,j);
%        plaxis(3,j)= 0;                %plaxis(4,j)
         end
      elseif iforward==5,
         plot(a(:,j1+1),a(:,j1+2),'b-',a(:,j1+3),a(:,j1+4),'r--')
      end
      xrev = 0;  yrev = 0;
      if (plaxis(4,j) < plaxis(3,j))
         xtemp = plaxis(:,j);
         plaxis(3,j) = xtemp(4);
         plaxis(4,j) = xtemp(3);
         yrev = 1;
      end
      if (plaxis(2,j)< plaxis(1,j))
         xtemp = plaxis(:,j);
         plaxis(1,j) = xtemp(2);
         plaxis(2,j) = xtemp(1);
         xrev = 1;
      end
      axis(plaxis(:,j));
      if (xrev==1); set(gca,'xdir','rev'); end
      if (yrev==1); set(gca,'ydir','rev'); end
      set(gca,'linewidth',1.2)
%%%
%  Put the xlabels, and ylabels for the comparisons of fields
%%%
      if rem(jplot,mhoriz) == 1 
         jtemp = 1 + jtemp;
         ylabels = char(title_str{jtemp},'Element #');
         ylabel(ylabels);
      end
      if jplot <= mhoriz
      elseif (jplot > mhoriz*(mvert-1)) 
         if rem(jplot,mhoriz) ~= 0
            xlabel(xlabel_str{rem(jplot,mhoriz)})
         else
            xlabel(xlabel_str{3})
         end
      end
%---------------------------------------------------
      xtext = ['curve ',num2str(j)];
      set(gca,'Fontsize',10);
   end
   if mvert >= 3 
      newpp = [1.25 2 6 7]; 
      set(gcf,'Paperposition',newpp) 
   end
end
%%%
%  Put comments and date on each figure:
%%%

if(exist('stringofcomments') == 0);
  stringofcomments0 = ['Results from data ',filenm,', (',num2str(bestfit),...
	' bestfit)'];
else
  stringofcomments0 =  stringofcomments;
end;

if(exist('no_title') == 0)
   figure(2)
   axes('position',[0,0,1,1]); axis('off');
   text(.05,.035, [stringofcomments0],'fontsize',9);
   text(.05,.02, ['File: ',pwd,'/',filenm,'; Date: ',date],'fontsize',6);
   filenmc = ['comparison',num2str(bestfit)];
   saveas(gcf,[filenmc '.fig'],'fig') 
   print('-depsc',[filenmc '.eps'])

   figure(1)
   axes('position',[0,0,1,1]); axis('off');
   text(.05,.035, [stringofcomments0],'fontsize',9);
   filenmp = ['ppd',num2str(bestfit)];
   text(.05,.02, ['File: ',pwd,'/',filenm,'; Date: ',date],'fontsize',6);
   saveas(gcf,[filenmp '.fig'],'fig') 
   print('-depsc',[filenmp '.eps'])  
end;


icom = true ;
%%%
% Data check: plotting the magnitudes of the data and the best model
%%%
if icom
   figure(3); clf
   jplot=0;
   for j = 1:nstep_plots:nplots
      jplot=jplot+1;
      if (rem(jplot,3)==0) 
         %   subplot(mvert,mhoriz,jplot)
         subplot(1,mvert,jplot/3)
         j1=4*(j-1);
         if (iforward==1 | iforward==2 | iforward==3|...
             (iforward==7&iopt(5)==4) | (iforward==8 & iopt(5)==4)|...
             iforward==9 ),
            plot(a(:,j1+1),a(:,j1+2),'b-',a(:,j1+3),a(:,j1+4),...
                 'r--','linewidth',1.2)
            bartarray_msr(:,jplot/3) = a(:,j1+1);
            bartarray_pre(:,jplot/3) = a(:,j1+3);
            if (plaxis(4,j)< plaxis(3,j));
               dum= plaxis(4,j);
               plaxis(4,j)= plaxis(3,j);
               plaxis(3,j)= dum;
               %  set(gca,'Ydir','reverse');
            end
            %      keyboard
         elseif ( (iforward==7 | iforward==8) &   iopt(5)==5), 
            plot(a(:,j1+1),a(:,j1+2),'b-',a(:,j1+3),a(:,j1+4),'r--')
         elseif (iforward==4| iforward==7 | iforward==8 ),
            plot(a(:,j1+1),a(:,j1+2),'b-',a(:,j1+3),a(:,j1+4),'r--')
         elseif iforward==6,    %TPEM
            if iopt(3)==12,
               plot(a(:,j1+2),a(:,j1+3),'b-',a(:,j1+4),a(:,j1+3),'r--')
               xtemp=plaxis(:,j);
               plaxis(1,j)=xtemp(3);
               plaxis(2,j)=xtemp(4);
               plaxis(3,j)=xtemp(1);
               plaxis(4,j)=xtemp(2);       
            else
               plot(a(:,j1+1),a(:,j1+2),'b-',a(:,j1+3),a(:,j1+4),'r--')
               dum= plaxis(4,j);
               plaxis(4,j)= plaxis(3,j);
               plaxis(3,j)= dum;
            end
         elseif iforward==5,
            plot(a(:,j1+1),a(:,j1+2),'b-',a(:,j1+3),a(:,j1+4),'r--')
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
         set(gca,'Fontsize',12);
      end
      set(gca,'linewidth',1.4,'Fontsize',12);
      if (rem(jplot,3)==0) 
         ylabel('Element #');
      end
      if (rem(jplot,3)==0)
         title(title_str{jplot/3})
      end
      set(gca,'ytick',[0:1:16])
      set(gca,'yticklabel','|1|2|3|5|6|7|8|9|10|11|12|13|14|15|16|')
      set(gca,'xlim',[0 0.3])
   end
   newpp = [1.25 3.5 7 4]; 
   set(gcf,'Paperposition',newpp) 
%   print -deps comparsion.eps
else
end