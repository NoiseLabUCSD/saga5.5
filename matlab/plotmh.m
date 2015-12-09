function  plotmh(filenm)
% use two Markov chains to do inference
ks = false;
%ks = true;
filename=[filenm '.mh1'];
read_mh_bin;
fitval=fval;
samp = xtt; 
if ks, samp1 = xtt; end
names = dir(filename);
% only read the *mh2 file if the mmatlab has enough room
if names.bytes < 90E6;
   filename = [filenm '.mh2'];
   read_mh_bin;
   fitval = [fitval; fval];samp=[samp; xtt];
 else
   display('Warning *mh2 file not read!')
   keyboard
 end
if ks, samp2 = xtt;end;
fval = []; xtt = [];
iplot = 2;
HPD = [50 75 95];
[bestofall bbin]= min(fitval);  % find the ML solution
aa = colormap(gray(64));
mycolormap = aa([1 24 36 50 64],:);
eval([filenm '_mh']);
iforward = iopt(30);
if (1==-1)
  figure
  % marginals from fortran program.
  for ii=1:nparm
    subplot(nparm ,1,ii)
    plot(x(:,ii),ppd(:,ii))
    %  set(gca,'xlim',[f_min(ii) f_max(ii)])
    xlabel( xtitles(iforward,par2phy(ii)),'Fontsize',10);
  end
end

npts=length(fitval)
maxl=f_max-f_min;
nbin=50;
bin=ones(nbin,1)*f_min'+[0.5:(nbin-0.5)]'*maxl'/(nbin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute parameter covariance matrix and correlation coefficient
% matrix
if (1==-1) 
   for i = 1:nparm
      mhist2(i,:) = (samp(:,i)'-f_min(i))/maxl(i);
      mhistm(i,:) = mhist2(i,:)-sum(mhist2(i,:))/npts;
   end
   for i = 1:nparm
      for j = 1:nparm
         cov(i,j) = sum(mhistm(i,1:npts).*mhistm(j,1:npts));   
      end
   end   
   cov=cov/npts;
   for i = 1:nparm
      for j = 1:nparm
         corr(i,j) = cov(i,j)/sqrt(cov(i,i)*cov(j,j));
      end
   end 
   % parameter covariance matrix
%   figure
%   waterfall(corr)
%   axis ij
%   view([0 78])
   save res1.mat corr f_min f_max cov
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = 1:nparm
   [temp1,temp2] = hist(samp(:,j),bin(:,j));
   binsamp(:,j) = temp1'/npts;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ks
   kk = 500; 
   aa = round(2.^linspace(7,floor(log2(size(samp1,1))),kk));
   ks = nan*ones(1,kk);chi2 = ks;
for ik = 1:kk
   for j = 2:2;%nparm
      if 1 == 1
        [H,P,ksstat] = kstest2(samp1([1:aa(ik)],j),...
                               samp2([1:aa(ik)],j));
        ks(j,ik) = ksstat;
        PP(j,ik) = P;
      else
      [temp11,temp2] = hist(samp1([1:aa(ik)],j),bin(:,j));
      bb = cumsum(temp11.');
      subsamp1 = bb./bb(end);
      [temp12,temp2] = hist(samp2([1:aa(ik)],j),bin(:,j));
      bb = cumsum(temp12.');
      subsamp2 = bb./bb(end);
      chi = sum(abs(temp12./aa(ik)-temp11./aa(ik)).^2);
      dif = abs(subsamp1-subsamp2);
      chi2(j,ik) = chi;
      ks(j,ik) = max(dif);
      end
   end
end
if size(ks,1) >= 10;
   n = 2;
else
   n = 1;
end
m = ceil(size(ks,1)/n);

figure;
orient landscape
for j = 1:size(ks,1)
   subplot(m,n,j)
   semilogx(aa,ks(j,:));hold on
   semilogx(aa,chi2(j,:),'r');
   semilogx([aa(1) aa(end)],[e_stop e_stop],'--r')
   set(gca,'ylim',[0 1])
   set(gca,'xlim',[aa(1) aa(end) ])
   text('position',[aa(end)/100*99 0.32],...
        'string', xtitles(iforward,par2phy(j)),...
        'horizontalalignment','right');
   if j == 1
      legend('KS','\chi^2')
   end
   if j == size(ks,1)
   xlabel('Number of samples')
   end
end 
  
figure(1)
   set(0,'DefaultAxesColorOrder',clrscl('rgb',13))
   semilogx(aa,ks.')
   ylabel('KS statistic')
   xlabel('Number of samples')
   style1(16,1.5)
   legend('sr','sd','wd','b','\theta','c','\Delta c','d','\alpha','\rho',...
       'EOF1','EOF2','EOF3')
end

ymax = max(max(binsamp));

if (1==-1)
  % timeseries
  figure
  for ii=1:nparm
    subplot(nparm,1,ii)
    plot(samp([1:round(end/2)],ii),'.','markersize',1)
    %  set(gca,'xlim',[f_min(ii) f_max(ii)])
    ylabel(xtitles(iforward,par2phy(ii)),'Fontsize',10);
    set(gca,'ylim',[f_min(ii) f_max(ii) ]); %axis ...
    axis tight
  end
end
%Index
nobs=npts;
ndigit=ones(nparm,1)*25; % the number of digits used in the plotting.
df=maxl./(ndigit);
parmind=uint8(round((samp(:,1:nparm)-ones(nobs,1)...
                     *f_min').*(ones(nobs,1)*(1./df'))+0.5));
% to avoid having any 0.
parmind(find(parmind<=0))=ones(size((find(parmind<=0))));
for ii=1:nparm
   parmind(find(parmind(:,ii)>ndigit(ii)),ii)=...
       ndigit(ii)* ones(size(find(parmind(:,ii)>ndigit(ii))));
   ml(ii) = samp(bbin,ii);
end

%save post.mat samp
 
% 2D marginals:
xval=zeros(nparm,max(ndigit));
for ii=1:nparm
  xval(ii,1:ndigit(ii))=([1:ndigit(ii)]-0.5)*df(ii)+f_min(ii);
end
% calulate 1-D marginals
marg=zeros(nparm,max(ndigit));
%keyboard

for ii=1:nparm
   for iobs=1:nobs
      ind= parmind(iobs,ii); 
      marg(ii,ind)= marg(ii,ind)+1;
   end
  [a b] = max(marg(ii,:));
  map(ii) = xval(ii,b);
end

marg=marg/nobs;
margmax=max(max(marg));

if (1==-1) 
   lcolor = 'b-';
   % MARGINALS
   figure(2);clf
   if par2phy(end)== 29; mparm = nparm-1; % i don't want to plot ppd of nu 
   else mparm = nparm
   end
   cols = 5;
   rows = ceil(mparm/cols);
   width = 0.9/cols;
   height = 0.9/rows;
   space = .2; % 2 percent space between axes
   for ii = 1:rows
      for jj = 1:cols
         kk =  (ii-1)*cols+jj;  
         if kk > mparm; return
         end
         axPos = [(jj-1)*width+0.5*space (rows-ii)*height+0.5*space ...
                  width*(1-space) height*(1-1.5*space)];
         subplot('position',axPos)
         plot(bin(:,kk), binsamp(:,kk),[lcolor],'linewidth',1);    
         hold on
         plot([ml(kk) ml(kk)],[-0.01 ymax*0.85],['-^' lcolor(1)],...
              'linewidth',1,'markerfacecolor', lcolor(1))
         set(gca,'linewidth',1.5)
         axis([f_min(kk) f_max(kk) 0 ymax*1.1]); %axis tight
         if exist('xtitle_user.m')
            xlabel(xtitle_user(kk),'Fontsize',10);      
         else
            xlabel(xtitles(iforward,par2phy(kk)),'Fontsize',10);
         end
         if jj > 1;
            set(gca,'yticklabel','');
         else
            ylabel('pdf','Fontsize',10);      
         end
      end      
      if mparm >= 4  
         newpp = [0.25 0.25 8 6]; 
         set(gcf,'Paperposition',newpp) 
      end
   end
   save marg1d.mat bin binsamp ml f_min f_max ymax par2phy iforward ...
       nparm mparm
end
save marg1d.mat bin binsamp ml f_min f_max par2phy iforward

% calulate 2-D  marginals
marg2d=zeros(nparm,nparm,max(ndigit),max(ndigit));
for ii=1:nparm
  for ii2=ii+1:nparm
    xx=zeros(max(ndigit),max(ndigit));
    for iobs=1:nobs
      ind= parmind(iobs,ii);
      ind2= parmind(iobs,ii2);
      xx(ind,ind2)= xx(ind,ind2)+1; %fitval(iobs);
    end
    marg2d(ii,ii2,:,:)= xx;
  end
end
marg2d=marg2d/nobs;
marg2dmax=max(max(max(max(marg2d))));

%2D filter
n1=3; sigma1=2; n2=3; sigma2=2; theta=0;
filter1 = d2gauss(n1,sigma1,n2,sigma2,theta);

if nparm > 10
   fontsize1d = 6;
   linew = 1;
else
   fontsize1d = 12;
   linew = 1.75;
end
figure('Name','1-D and 2-D PPDs',...%'NumberTitle','off',...
       'PaperPosition',[0.25 2.5 8 6]);
pos = [0.05 0.1 0.875 0.85];
axes('position',pos); axis('off')

%keyboard
space = .1; % 10 percent space between axes
%width = (pos(3)-space)/nparm;
%height = (pos(4)-2.5*space)/nparm;
width = (pos(3))/nparm;
height = (pos(4))/nparm;
pos(1:2) = pos(1:2) + space*[width height];

% 2D marginals
nparmm=nparm-1;
iflag = 0;
for ii=1:nparmm
   for ii2=ii+1:nparm
      axPos = [(ii2-1)*width+pos(1) (nparm-ii)*height+pos(2) ...
               width*(1-space) height*(1-space)];
      subplot('position',axPos);
      zcontr = squeeze(marg2d(ii,ii2,1:ndigit(ii),1:ndigit(ii2))); 
      zcontr = conv2(zcontr,filter1,'same')*100;
      zcontr1 = zcontr;
      if iplot == 1
         imagesc(xval(ii2,1:ndigit(ii2)),xval(ii,1:ndigit(ii)),...
                 zcontr,[0 marg2dmax]);axis xy
         axis([f_min(ii2) f_max(ii2) f_min(ii) f_max(ii)]);
         if iflag == 0;
            iflag = 1;
            lhand = colorbar('vert','XAxisLocation','top',...
                             'fontsize',14);
            xlabel(lhand,'pdf')
         end
      else 
         if isempty(HPD)
            ZHPD   = 5; % contour lines
            contourf(xval(ii2,1:ndigit(ii2)),...
                     xval(ii,1:ndigit(ii)),...
                     zcontr,ZHPD,'edgecolor','none');
            %caxis([0 0.4*marg2dmax]);%%%%%%%%%%%%%%%%
            %caxis([0 0.8*marg2dmax]);%%%%%%%%%%%%%%%%
            colormap(flipud(pink(32)))
            if iflag == 0;
               iflag = 1;
               lhand = colorbar('vert','XAxisLocation','top',...
                                'fontsize',14);
               xlabel(lhand,'pdf')
            end
         else
            ZHPD  = findcontour(zcontr,HPD);
            zhpdma= max(ZHPD);
            [c,h,cf] = contourcf(xval(ii2,1:ndigit(ii2)),...
                                xval(ii,1:ndigit(ii)),...
                                zcontr./zhpdma,[ZHPD 0]./zhpdma);
            %set(h,'edgecolor','none') 
            colormap(flipud(mycolormap))
            oldh = get(h,'FaceVertexCData');  
            if iflag == 1
            else
               oldh1 = cell2mat(oldh);
               [ia ib] = sort(oldh1); 
               if isempty(find(diff(ia) == 0)) & (iflag == 11);
                  iflag = 1;
                  lhand = legend(...
                      [h(end) h(end-1) h(end-2)],...
                      [num2str(HPD(1)) '% HPD'],...
                      [num2str(HPD(2)) '% HPD'],...
                      [num2str(HPD(3)) '% HPD']);
               end
            end
         end
         [zi,zj] = find(zcontr1 == max(max(zcontr1)));
         hold on
         if nparm > 10
            plot(xval(ii2,zj), xval(ii,zi),'+w','markersize',6,...
                 'linewidth',1); % mode
            plot(ml(ii2),ml(ii),'wx','markersize',6,'linewidth',1);
         else
            plot(xval(ii2,zj), xval(ii,zi),'+w','markersize',8,...
                 'linewidth',1.5); % mode
            plot(ml(ii2),ml(ii),'wx','markersize',8,'linewidth',1.5);
         end
         axis([f_min(ii2) f_max(ii2) f_min(ii) f_max(ii)]);axis xy
      end
      if (rem((ii-1)*nparm+ii2,nparm) == 0 &...
         ((ii-1)*nparm+ii2)/nparm >= 1)
         set(gca,'XTickLabel','','YAxisLocation','Right')
      else
         set(gca,'XTickLabel','','YTickLabel',' ')
      end
      if nparm > 5
         set(gca,'XTickLabel','','YTickLabel',' ')
      end
      style(fontsize1d,linew);
   end
end

for ii=1:nparm
   axPos = [(ii-1)*width+pos(1) (nparm-ii)*height+pos(2) ...
            width*(1-space) height*(1-space)];
   subplot('position',axPos);
   plot(xval(ii,1:ndigit(ii)),marg(ii,1:ndigit(ii)),'linewidth',1)
   hold on
   plot([ml(ii) ml(ii)],[-0.01 margmax*0.9],'-^r',...
        'linewidth',1,'markerfacecolor','r')
   if exist('truemodel.m')
      if ~isempty(truemodel(ii))
         plot([truemodel(ii) truemodel(ii)],[0 margmax],':k',...
           'linewidth',1)
      end
   end
   style(fontsize1d,linew)
   axis([f_min(ii) f_max(ii) 0 margmax])
   if exist('xtitle_user.m')
      xlabel(xtitle_user(ii),'Fontsize', fontsize1d);      
   else
      xlabel(xtitles(iforward,par2phy(ii)),'Fontsize',fontsize1d);
   end
   if ii > 1
      set(gca,'YTickLabel',' ')
   end
end

if ~isempty(HPD)
  % set(lhand,'position',[0.08 0.45 0.15 0.0975])
else      
   set(lhand,'position',[0.13 0.275 0.025 0.15])
end
axes('position',[0,0,1,1]); axis('off')
xx = text(.05,.02, ['File: ',pwd,'/',filenm,'; Date: ',date],...
          'fontsize',6,'interpreter','none');

%     
%---- Estimating standard deviation
%     
for ii=1:nparm
   xsum=0.;
   xsumsqr=0.;
   for jj=1:ndigit(ii)
      xhelp = xval(ii,jj).*marg(ii,jj);
      xsum = xsum+xhelp;
      xsumsqr = xsumsqr+xhelp.*xval(ii,jj);
   end
   xmean(ii)=xsum;
   sstd(ii)=sqrt(abs(xsumsqr-xsum.^2));
end
save statist.mat xmean sstd
%save res.mat fitval

if (1==-1)
   %convergence
   figure
   eval(['load -ascii  ' filenm '.conv']) 
   eval(['convdata=' filenm '(:,3:end) ;']);
   for ii=1:nparm
      subplot(nparm,1,ii)
      semilogx([1:size(convdata,1)], convdata(:,ii),'-')
      hold on
      semilogx([1 size(convdata,1)],[e_stop e_stop],'--')
      axis([1 size(convdata,1) 0 1])
      text('position',[(length(convdata(:,1)))/100*99 0.32],...
           'string', xtitles(iforward,par2phy(ii)),...
           'horizontalalignment','right');
   end
   xlabel('# of runs')
end  

