function  plotgibbs(filenm)

%  keyboard
filename=[filenm '.gs1'];
read_gs_bin;
fitval=fval;
samp=xtt;
filename=[filenm '.gs2'];
read_gs_bin;
fitval=[fitval; fval];
samp=[samp; xtt];



eval([filenm '_gib']);
iforward= iopt(30);
if (1==-1)
  figure
  % marginals from fortran program.
  for ii=1:nparm
    subplot(6,1,ii)
    plot(x(:,ii),ppd(:,ii))
    %  set(gca,'xlim',[f_min(ii) f_max(ii)])
    xlabel( xtitles(iforward,par2phy(ii)),'Fontsize',10);
  end
end

%samp=[samp1; samp2];

%data=samp(:,2:(nparm+1));
%fit=samp(:,nparm+2);
npts=length(fitval)
maxl=f_max-f_min;
nbin=50;
bin=ones(nbin,1)*f_min'+[0.5:(nbin-0.5)]'*maxl'/(nbin);
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
for j = 1:nparm
   [temp1,temp2] = hist(samp(:,j),bin(:,j));
   binsamp(:,j) = temp1'/npts;
end

ymax=max(max(binsamp));

if (1==-1) 
% MARGINALS
  figure
  for ii=1:nparm
    subplot(6,1,ii)
    plot(bin(:,ii), binsamp(:,ii))
    %  set(gca,'xlim',[f_min(ii) f_max(ii)])
    xlabel( xtitles(iforward,par2phy(ii)),'Fontsize',10);
    axis([f_min(ii) f_max(ii) 0 ymax*1.1]); %axis tight
  end
end

if (1==-1)
  % timeseries
  figure
  for ii=1:nparm
    subplot(6,1,ii)
    plot(samp(:,ii),'.','markersize',1)
    %  set(gca,'xlim',[f_min(ii) f_max(ii)])
    ylabel( xtitles(iforward,par2phy(ii)),'Fontsize',10);
    set(gca,'ylim',[f_min(ii) f_max(ii) ]); %axis ...
    tight
  end
end
%Index
nobs=npts;
ndigit=ones(nparm,1)*30; % the number of digits used in teh plotting.
df=maxl./(ndigit);
parmind=uint8(round((samp(:,1:nparm)-ones(nobs,1)*f_min').*(ones(nobs,1)*(1./df'))+0.5));
% to avoid having any 0.
parmind(find(parmind<=0))=ones(size((find(parmind<=0))));
for ii=1:nparm
parmind(find(parmind(:,ii)>ndigit(ii)),ii)=...
       ndigit(ii)* ones(size(find(parmind(:,ii)>ndigit(ii))));
end
% 2D marginals:
xval=zeros(nparm,max(ndigit));
for ii=1:nparm
  xval(ii,1:ndigit(ii))=([1:ndigit(ii)]-0.5)*df(ii)+f_min(ii);
end
% calulate 1 D  marginals
marg=zeros(nparm,max(ndigit));
%keyboard
for ii=1:nparm
for iobs=1:nobs
    ind= parmind(iobs,ii); %fix((xtt(iobs,ii)-f_min(ii))/df(ii))+1;
    marg(ii,ind)= marg(ii,ind)+1; %fitval(iobs);
end
end
marg=marg/nobs;
margmax=max(max(marg));


% calulate 2-D  marginals
marg2d=zeros(nparm,nparm,max(ndigit),max(ndigit));
for ii=1:nparm
  for ii2=ii+1:nparm
    xx=zeros(max(ndigit),max(ndigit));
    for iobs=1:nobs
      ind= parmind(iobs,ii); %fix((xtt(iobs,ii)-f_min(ii))/df(ii))+1;
      ind2= parmind(iobs,ii2); % fix((xtt(iobs,ii2)-f_min(ii2))/df(ii2))+1;
      xx(ind,ind2)= xx(ind,ind2)+1; %fitval(iobs);
    end
    marg2d(ii,ii2,:,:)= xx;
  end
end
marg2d=marg2d/nobs;
marg2dmax=max(max(max(max(marg2d))));

%

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

if (1==-1)
%convergence
figure
eval(['load -ascii  ' filenm '.conv']) %samp2.mat   %[filenm '.gs1']
eval(['convdata=' filenm '(:,3:end) ;']);

%keyboard
for ii=1:nparm
  subplot(6,1,ii)
 %  plot(convdata(:,ii),'-')
   semilogy(convdata(:,ii),'-')
hold on
semilogy([0 length(convdata(:,1))],[ 0.05 0.05],'--')
 axis([0 length(convdata(:,1)) 0.03 0.5])
end
end  
  %keyboard