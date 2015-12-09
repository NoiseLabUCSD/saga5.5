function [xval,marg]= plotenum(filenm,indx)
% plotting marginal distribution from the SAGA enumerative integration
% module
if nargin==1,
  indx=10
end

  filename=[filenm '.enum'];
%filename='enum.mat'
read_mh_bin;
fitval = fval; % Note: fval = func/nugibbs
[bestofall bbin]= min(fval);
samp=xtt;
iforward=iopt(30);
npts = length(fitval)
maxl=f_max-f_min;

%Index
nobs = npts;
parmind=uint8(round((samp(:,1:nparm)-ones(nobs,1)*f_min').*...
                    (ones(nobs,1)*(1./df')))+1);
% to avoid having any 0.
parmind(find(parmind<=0)) = ones(size((find(parmind<=0))));
for ii=1:nparm
   parmind(find(parmind(:,ii)>ndigit(ii)),ii) = ndigit(ii)* ...
       ones(size(find(parmind(:,ii)>ndigit(ii))));
end

expfit = exp(-fitval);
mdigit = max(ndigit);
%2D marginals:
xval=zeros(nparm,max(ndigit));
for ii=1:nparm
   xval(ii,1:ndigit(ii))=([1:ndigit(ii)]-1.0)*df(ii)+f_min(ii);
end

% calulate 1D marginals
marg=zeros(nparm,max(ndigit));
for ii=1:nparm
   fac = ndigit(ii)/mdigit;
   for iobs=1:nobs
      ind= parmind(iobs,ii); %fix((xtt(iobs,ii)-f_min(ii))/df(ii))+1;
      marg(ii,ind)= marg(ii,ind)+ fac*expfit(iobs);
   end
%   [a b] = max(marg(ii,:));
   map(ii) = xtt(bbin,ii);
end
marg = marg/sum(expfit);
margmax = max(max(marg));

% calulate 2D marginals
marg2d = zeros(nparm,nparm,max(ndigit),max(ndigit));
for ii = 1:nparm
   fac1 = ndigit(ii)/mdigit;
   for ii2 = ii+1:nparm
      fac2 = ndigit(ii2)/mdigit;
      xx = zeros(max(ndigit),max(ndigit));
      for iobs=1:nobs
         ind  = parmind(iobs,ii); %fix((xtt(iobs,ii)-f_min(ii))/df(ii))+1;
         ind2 = parmind(iobs,ii2); % fix((xtt(iobs,ii2)-f_min(ii2))/df(ii2))+1;
         xx(ind,ind2) = xx(ind,ind2) + fac1*fac2*expfit(iobs);
      end
      marg2d(ii,ii2,:,:) = xx;
   end
end
marg2d = marg2d/sum(expfit);
marg2dmax = max(max(max(max(marg2d))));

% keyboard;
thisfig=gcbf;
if ishandle(thisfig)
   nextfig=thisfig+1;
else
   nextfig=1;
end

figure(indx);
% 2D marginals
nparmm=nparm-1;
kk = 0;
for ii=1:nparmm
   for ii2=ii+1:nparm
      kk = kk+1;
      ax(kk) = subplot(nparm,nparm,(ii-1)*nparm+ii2);
      imagesc(xval(ii2,1:ndigit(ii2)),xval(ii,1:ndigit(ii)),...
              squeeze(marg2d(ii,ii2,1:ndigit(ii),1:ndigit(ii2))),...
              [0 0.7*marg2dmax]);
      axis xy
   end 
end
if (nparm>1)
  hh = colorbar;
  pos = get(ax(nparmm), 'Position');
  poskk = get(ax(kk), 'Position');
  set(ax(kk), 'Position', [poskk(1) poskk(2) pos(3) poskk(4)])    
  set(hh,'Position',...
      [(pos(1)+pos(3))*1.02 pos(2) .01 pos(4)*0.8])
end
  
for ii=1:nparm
   subplot(nparm,nparm,(ii-1)*nparm+ii)
   plot(xval(ii,1:ndigit(ii)),marg(ii,1:ndigit(ii)),'b','linewidth',1);
   hold on
   set(gca,'linewidth',1)
   plot([map(ii) map(ii)],[-0.01 margmax*0.9],'-^r',...
        'linewidth',1,'markerfacecolor','r')
   set(gca,'xlim',[f_min(ii) f_max(ii)],'ylim',[0 margmax])
   xlabel(xtitles(iforward,par2phy(ii)),'Fontsize',10);
   map(ii)
 end

% for standard deviation
for ii=1:nparm
   xsum=0.;
   xsumsqr=0.;
   for jj=1:ndigit(ii)
      xhelp = xval(ii,jj).*marg(ii,jj);
      xsum = xsum+xhelp;
      xsumsqr = xsumsqr+xhelp.*xval(ii,jj);
   end
   xmean(ii)=xsum;
   svar(ii)=sqrt(abs(xsumsqr-xsum.^2));
end


eval(['save res', num2str(indx)])

%keyboard
