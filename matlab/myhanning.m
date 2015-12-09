function [han]=hanning(fmin,fmax,freq,df);
%
%   This function produces a Hanning window
%
fc=(fmax+fmin)/2;
%
han=zeros(1,size(freq,2));
ind1=ceil(fmin/df)+1;
ind2=fix(fmax/df)+1;
han(ind1:ind2)=cos(pi*(freq(ind1:ind2)-fc)/(fmax-fmin));
%
%plot(freq,han)
