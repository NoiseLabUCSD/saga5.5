% COSTAP(n,nf1,nf2,percent)  Cosine taper window 
% The window can be placed anywhere within the sequence
% use nf1=1, nf2=n for full-length cosine taper window
%       n - total window length
%       nf1 - lower index of tapered window (start taper)
%       nf2 - upper index of tapered window (stop taper)
%       percent - taper length at each end as a percent of n  

% in order of increasing magnitude: 1, nf1, n1, n2, nf2, n.
% 1-f1: 0's
% nf1-n1: up taper
% n1-n2: 1's
% n2-nf2: down taper
% nf2-n: 0's

function w = costap(n,nf1,nf2,percent)


ntap=floor((nf2-nf1+1)*percent/100); % number of points in taper
n1=nf1+ntap; % first point after tapered section ends
n2=nf2-ntap; % last point before down taper section starts
dpi=3.141592654 /(ntap-1); % pi over ntap-1 intervals in taper

w=ones(1,n); % allocate the memory and initialize all to 1
w(1:nf1)= zeros(1,nf1); 
w(nf1:n1-1)=(1-cos([0:ntap-1]*dpi ))*0.5;
w(nf2:-1:n2+1) = w(nf1:n1-1);
w(nf2:n)=zeros(1,n-nf2+1); 
                
return


