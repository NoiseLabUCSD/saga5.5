
for idata = 1:ndata
idata
icsdm = datablocks(idata);
   eval(['load ',datapath,'csdm',num2str(nsnap),'_',num2str(icsdm)]);
   disp(['Time = ',num2str(Time(icsdm)/60,3),' min']); 
   [nele, ncov] = size(Rf);
   nfre        = ncov/nele;
   Cov         = reshape(Rf, nele, nele, nfre);
   K           = Cov(:,:,ifre);
   ktrace      = trace(K)
   K           = K ./ trace(K); % Finally, we have data csdm
   
   %--- estimate snr
%   [EIG Dig V] = svd(K);  
%   D           = diag(Dig);
%   snr         = D(1)/nele;
   
   
   % initialization
   %--- compute the likelihood
% depth search band
 isband=2;
 idb=idepth+[-isband:isband];
 ii=find(idb<1);  idb(ii)=ones(size(ii)); 
     ii=find(idb>20);  idb(ii)=20*ones(size(ii)); 
 ienv2=ones(2*isband+1,1)*[1:nenv]+(idb-1)'* nenv*ones(1,nenv);
ienv=reshape((ienv2)',nenv*(2*isband+1),1);
% search band for sd
%%% search band for range
iband=5;
    irangeint=irange+ [-iband:iband];
     ii=find(irangeint<1);  irangeint(ii)=ones(size(ii)); 
     ii=find(irangeint>par(7));  irangeint(ii)=par(7)*ones(size(ii)); 

   for iii = 1:length(irangeint)
     irng =irangeint(iii);
     repl1 = squeeze(s1(:,irng,ienv));
      for ii = 1:length(ienv) %1:menv
	 ienv1=ienv(ii);
         repl =  repl1(:,ii);
%        repl =  repl1(:,ienv1);
         BB(irng,ienv1) = abs(repl' * K * repl) ;%/norm(repl)^2;
      end
    end
%%% search band for range
[dd,ii] = max(BB(irangeint,ienv),[],2);  %Bmax(idata) = max(max(BB(irangeint,:)));         % find ML solution 
 [dd,ir] = max(dd);
maxenv= ienv(ii(ir));

 irange=irange+ir-iband-1;
 Bmax(idata)=dd;
maxrng=irange
 % [maxrng, maxenv] = find(BB == Bmax(idata)); 
   sr(idata) = ranges(maxrng);              % ML sr
   isd = floor(maxenv/(nenv));
   idepth=isd+1
   sd(idata)  = deptrial(isd+1);            % ML sd
   disp(['MEASURED (sr, sd) = (',num2str(sr_save(icsdm)),',',...
         num2str(sd_save(icsdm)),')']); 
   disp(['ESTIMATED (sr, sd) = (',num2str(sr(idata)),',', num2str(sd(idata)),')']); 
   nu = (1-Bmax(idata))/Neff;
   obj(idata,:) = (nu^(-Neff))*...
       exp(-(1-abs(BB(maxrng,[nenv*isd+1:nenv*(isd+1)])))./nu);
   phi(idata,:)=ktrace*(1-abs(BB(maxrng,[nenv*isd+1:nenv*(isd+1)])));
 end
