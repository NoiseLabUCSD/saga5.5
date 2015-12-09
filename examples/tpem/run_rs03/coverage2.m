function [OneFreq] = coverage2(s,par,ifreq,crange,xplot);
%COVERAGE	picks the ifreq frequency from TPEM data,
%               and plots the coverage diagram for the particular frequency.
%		First read it using READASCII
%
%	Syntax:
%
%	Onefreq=coverage(s,par,ifreq;
%
%	Onefreq	-  complex frequency response in range and depth
%	
%	s,par	- input variable from READASCII
%       ifreq	- frequency to plot
%       crange  - range for color plot

%	receiv	- vector with wanted receivers
%		  (scalar for only one)
%	range	- vector with wanted ranges
%		  (scalar for only one)


%----------------
nout	= par(1);
rd	= par(2);
rdlow	= par(3);
ir	= par(4);
r0	= par(5);
rspace	= par(6);
nplots	= par(7);
nx	= par(8);
lx	= par(9);
mx	= par(10);
dt	= par(11);
msuft	= par(12);
df=1/(par(11)*par(8));            %  delta freq
fmin=df*(par(9)-1);               % the minimum used frequency
fmax=df*(par(10)-1);              % the maximum used frequency
freq=df*((par(9)-1):(par(10)-1)); % all the used frequencies
%----------------
range   = 1:par(7);
receiv  = 1:par(4);
if (max(receiv)>ir)|(min(receiv)<1),
	error('receivers out of range');
end;

if (max(range)>nplots)|(min(range)<1),
	error('range out of range');
end;


dx	= msuft*nplots*ir;	
%f	= 0:dx:(mx-lx)*dx;
f	= (ifreq-1)*dx;
dmsuft	= nplots*ir;		
dnplots	= ir;			

sr	= [];
%keyboard
for x = (range-0),
	for r=receiv,
%		sr(par(4)+1-r,x)	=  s(f+(x-1)*dnplots+r,1);
		sr(r,x)	=  s(ifreq,(x-1)*dnplots+r);
%		sr(r,x)	=  s(f+(x-1)*dnplots+r,1);
	end;
end;

 height=linspace(rd,rdlow,ir);
 ranges=par(5)+par(6)*( 0:par(7)-1);
OneFreq=sr;
index= find(abs(sr)==0); 
%sr(index)=ones(size(index));  
%keyboard 
sr=20*log10(abs(sr));            %-ones(size(sr,1),1).*10*log10(range*1000);
sr(index)=1e9*ones(size(index));   
%pcolor(ranges,height,sr),shading('flat'),h=colorbar('vert'); 
axis('xy')
if (nargin==5)
  pcolor(ranges,height,sr), caxis(crange),shading('flat'),h=colorbar('vert');
elseif (nargin==4)
  imagesc(ranges,height,sr,crange),shading('flat'),h=colorbar('vert'); 
else
  imagesc(ranges,height,sr),shading('flat'),h=colorbar('vert'); 
end

%keyboard
% shading('interp')colorbar,colormap('gray'); 
pos=get(h,'position'); pos(3)=0.025; pos(1)=pos(1)-0.023;
set(h,'position',pos,'xlim',[0.3 1])
%
xlabel('Range (km)'),ylabel('Height (m)');
title(['Coverage diagram, freq ' int2str(freq(ifreq))])





