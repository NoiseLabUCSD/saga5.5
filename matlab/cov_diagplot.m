function cov_diagplot(freqbin,frequency,Cxest)
%cov_diagplot(freqbin,frequency,Cxest)
%  
%  pretty-print plot of characteristics of a covariance matrix
%  used by for debugging the signal processing
%  Christoph Meclenbrauker & Peter Gerstoft
subplot(2,3,1); 
     mesh(abs(Cxest));  title('abs(Cxest)')
subplot(2,2,2);
     plot(10*log10(svd(Cxest)),'o'); 
     title('SVD'); ylabel('Power [dB]');
subplot(2,3,4);
     imagesc(abs(Cxest)); xlabel('abs(Cxest)');
     s = sprintf('f = %10.4f Hz',frequency);
     title(s);
subplot(2,3,5); 
     imagesc(angle(Cxest)); xlabel('angle(Cxest)');
     s = sprintf('PG: freqbin = %d',freqbin);
     title(s);
return;


