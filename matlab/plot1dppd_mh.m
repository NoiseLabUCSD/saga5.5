function  plot1dppd_mh(filenm)

% function to plot the 1dppd of the parameters after running Markov Chain
% Monte Carlo in SAGA. filenm is the filename of the .gs1 file output from
% SAGA. f195hla_oS1gm

filename=[filenm '.mh1'];
names = dir(filename);
read_mh_bin;
fitval = fval;
samp = xtt; 
fval = []; xtt = [];

[bestofall bbin]= min(fitval);  % find the ML solution
eval([filenm '_mh']);
iforward = iopt(30);
npts = length(fitval);
nobs = npts;
maxl = f_max-f_min;
ndigit = ones(nparm,1)*25; % the number of digits used in the plotting.
df = maxl./(ndigit);
parmind = uint8(round((samp(:,1:nparm)-ones(nobs,1)...
                       *f_min').*(ones(nobs,1)*(1./df'))+0.5));
% to avoid having any 0.
parmind(find(parmind<=0)) = ones(size((find(parmind<=0))));
for ii=1:nparm
   parmind(find(parmind(:,ii)>ndigit(ii)),ii) =...
       ndigit(ii)* ones(size(find(parmind(:,ii)>ndigit(ii))));
   ml(ii) = samp(bbin,ii);
end

% 2D marginals:
xval = zeros(nparm,max(ndigit));
for ii = 1:nparm
  xval(ii,1:ndigit(ii))=([1:ndigit(ii)]-0.5)*df(ii)+f_min(ii);
end
% calulate 1-D marginals
marg = zeros(nparm,max(ndigit));
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


% Plot 1D-ppd
for ii=1:nparm
    figure;
    plot(xval(ii,1:ndigit(ii)),marg(ii,1:ndigit(ii)),'linewidth',1)
    hold on
    plot([ml(ii) ml(ii)],[-0.01 margmax*0.9],'-^r',...
        'linewidth',1,'markerfacecolor','r')
    axis([f_min(ii) f_max(ii) 0 margmax])
    xlabel( xtitles(iforward,par2phy(ii)),'Fontsize',18);    
    ylabel('pdf','Fontsize',18);
end

