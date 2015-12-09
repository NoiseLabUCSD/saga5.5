filename='enum.mat'
read_enum_bin

disp('... data read in')
xval=zeros(nparm,max(ndigit));
for ii=1:nparm
  xval(ii,1:ndigit(ii))=([1:ndigit(ii)]-1)*df(ii)+f_min(ii);
end

parmind=uint8(round((xtt(:,1:nparm)-ones(nobs,1)*f_min').*(ones(nobs,1)*(1./df')))+1);

% calulate 1 D  marginals
marg=zeros(nparm,max(ndigit));
for ii=1:nparm
for iobs=1:nobs
    ind= parmind(iobs,ii); %fix((xtt(iobs,ii)-f_min(ii))/df(ii))+1;
    marg(ii,ind)= marg(ii,ind)+fval(iobs);
end
end
margmax=max(max(marg));

disp('calculating 2-D  marginals')
% calulate 2-D  marginals
marg2d=zeros(nparm,nparm,max(ndigit),max(ndigit));
for ii=1:nparm
  for ii2=ii+1:nparm
    xx=zeros(max(ndigit),max(ndigit));
    for iobs=1:nobs
      ind= parmind(iobs,ii); %fix((xtt(iobs,ii)-f_min(ii))/df(ii))+1;
      ind2= parmind(iobs,ii2); % fix((xtt(iobs,ii2)-f_min(ii2))/df(ii2))+1;
      xx(ind,ind2)= xx(ind,ind2)+fval(iobs);
    end
    marg2d(ii,ii2,:,:)= xx;
  end
end
marg2dmax=max(max(max(max(marg2d))))
%
figure
for ii=1:nparm
  subplot(6,1,ii)
  plot(xval(ii,1:ndigit(ii)),marg(ii,1:ndigit(ii)))
  set(gca,'xlim',[f_min(ii) f_max(ii)],'ylim',[0 margmax])
end

figure
% 2D marginals
nparmm=nparm-1;

for ii=1:nparmm
  for ii2=ii+1:nparm
    subplot(nparm,nparm,(ii-1)*nparm+ii2)
    imagesc(xval(ii2,1:ndigit(ii2)),xval(ii,1:ndigit(ii)),...
    squeeze(marg2d(ii,ii2,1:ndigit(ii),1:ndigit(ii2))),[0 marg2dmax]),colorbar
  end
end
for ii=1:nparm
  subplot(nparm,nparm,(ii-1)*nparm+ii)
  plot(xval(ii,1:ndigit(ii)),marg(ii,1:ndigit(ii)))
  set(gca,'xlim',[f_min(ii) f_max(ii)],'ylim',[0 margmax])
%  xlabel( xtitles(iforward,par2phy(ii)),'Fontsize',10);
end



