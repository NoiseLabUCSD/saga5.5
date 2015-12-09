function [VCov,Signal,Noise] = svdanaly2(freq,Cov,Nshot,frequency)
%[VCov,Signal,Noise] = svdanaly(freq,Cov,Nshot,frequency)
%
%  Make a plot of the estimated SNR based on the 
%  SVD of the covariance over all frequencies.
%
[Nsensor,jmaxNsensor] = size(Cov);
jmax = jmaxNsensor/Nsensor;

fprintf(1,'svd SNR analysis\n');
VCov = zeros(size(Cov));

for jj=1:jmax
  jjvec    = [(Nsensor*(jj-1)+1):(Nsensor*jj)];
  Cxest    = Cov(:,jjvec);
  [U,S,V]  = svd(Cxest);
  VCov(:,jjvec) = V;
  Noise(jj)  = (trace(Cxest)-S(1,1))/(Nsensor-1);
  Signal(jj) = S(1,1);
%  Signal(jj) = S(1,1)-Noise(jj);
  if(1==0) 
    S = [S(1,1) Noise(jj) * ones(1,Nsensor-1)];
    Cov(:,jjvec) = U*diag(S)*V';
  end
end   

yvalue = 10*log10(Signal./Noise);
%yvalue = 10*log10(Signal)+13.98+138+20;
xvalue = freq;
coeff  = polyfit(xvalue,yvalue,2);

%xvalue1=linspace(min(freq),max(freq),100)';
xvalue1 = linspace(0.97*min(freq),1.03*max(freq),100)';

lsfit1 = polyval(coeff,xvalue1);

w=find(freq==frequency);

subplot(2,3,1),
plot([xvalue(1:w-1) xvalue(w+1:length(xvalue))],[yvalue(1:w-1) yvalue(w+1:length(yvalue))],'mo',xvalue1,lsfit1,'m',xvalue(w),yvalue(w),'bx');
set(gca,'fontsize',10);
title('a) Power','fontsize',12); 
xlabel('Frequency (Hz)','fontsize',10);
ylabel('Power (dB)','fontsize',10);
%grid

drawnow;

aux = lsfit1 - max(lsfit1);

MaximumIndex = find(aux==0);

fprintf(1,'Center Frequency = %f Hz\n',xvalue1(MaximumIndex));

ThreeDB = find(aux >= -3);

fprintf(1,'Lower 3dB cut-off = %f Hz\n',xvalue1(min(ThreeDB)));
fprintf(1,'Upper 3dB cut-off = %f Hz\n',xvalue1(max(ThreeDB)));

end;