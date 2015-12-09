function  plotenum(filenm)
% plotting marginal distribution from the SAGA enumerative integration module
filename=[filenm '.enum'];
%filename='enum.mat'
read_gs_bin;
fitval=fval;
samp=xtt;

iforward=iopt(30);
npts=length(fitval)
maxl=f_max-f_min;



%Index
nobs=npts
parmind=uint8(round((samp(:,1:nparm)-ones(nobs,1)*f_min').*(ones(nobs,1)*(1./df')))+1);
% to avoid having any 0.
parmind(find(parmind<=0))=ones(size((find(parmind<=0))));
for ii=1:nparm
  parmind(find(parmind(:,ii)>ndigit(ii)),ii)=...
       ndigit(ii)* ones(size(find(parmind(:,ii)>ndigit(ii))));
 end
 
expfit=exp(-fitval);
mdigit= max(ndigit);
% 2D marginals:
xval=zeros(nparm,max(ndigit));
for ii=1:nparm
  xval(ii,1:ndigit(ii))=([1:ndigit(ii)]-1.0)*df(ii)+f_min(ii);
end
% calulate 1 D  marginals
marg=zeros(nparm,max(ndigit));
for ii=1:nparm
for iobs=1:nobs
    ind= parmind(iobs,ii); %fix((xtt(iobs,ii)-f_min(ii))/df(ii))+1;
    marg(ii,ind)= marg(ii,ind)+expfit(iobs);
  end
end
marg=marg/sum(expfit);
margmax=max(max(marg));

%keyboard

% calulate 2-D  marginals
marg2d=zeros(nparm,nparm,max(ndigit),max(ndigit));
for ii=1:nparm
  for ii2=ii+1:nparm
    xx=zeros(max(ndigit),max(ndigit));
    for iobs=1:nobs
      ind= parmind(iobs,ii); %fix((xtt(iobs,ii)-f_min(ii))/df(ii))+1;
      ind2= parmind(iobs,ii2); % fix((xtt(iobs,ii2)-f_min(ii2))/df(ii2))+1;
      xx(ind,ind2)= xx(ind,ind2)+expfit(iobs);
    end
    marg2d(ii,ii2,:,:)= xx;
  end
end
marg2d=marg2d/sum(expfit);
marg2dmax=max(max(max(max(marg2d))));

%
%keyboard;

figure
% 2D marginals
nparmm=nparm-1;

for ii=1:nparmm
  for ii2=ii+1:nparm
    subplot(nparm,nparm,(ii-1)*nparm+ii2)
    imagesc(xval(ii2,1:ndigit(ii2)),xval(ii,1:ndigit(ii)),...
    squeeze(marg2d(ii,ii2,1:ndigit(ii),1:ndigit(ii2))),[0 0.7*marg2dmax] ),colorbar
    axis xy
end
end
for ii=1:nparm
  subplot(nparm,nparm,(ii-1)*nparm+ii)
  plot(xval(ii,1:ndigit(ii)),marg(ii,1:ndigit(ii)))
  set(gca,'xlim',[f_min(ii) f_max(ii)],'ylim',[0 margmax])
  xlabel( xtitles(iforward,par2phy(ii)),'Fontsize',10);
end


%keyboard
