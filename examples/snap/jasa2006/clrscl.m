function [ cmap]=clrscl(pal,ncol)
% clrscl.m
%-------------------
% The function clrscl is written by MAD and produces color scales
% that linear interpolote ncol colors between the points in color
% space. Currently it recognizes a colorwheel
% b blue
% c cyan
% g green
% y yellow
% r red
% m magenta
%
% also
% w white
% k black
%-------------------
% pal chooses route around color wheel.
% pal = 'wk' % white to black greyscale
% pal = 'rb' % red to blue
% pal = 'rbg' % red to blue to green
% pal = 'bcgyrm' % is a color wheel
% pal = 'rmbcgyr' % is a color wheel that is periodic, for phase.
% magenta (m) is between red and blue, cyan (c) is between blue and green,
% yellow (y) is between green and red.
%-------------------
% example:
% (to see the scale, don't give an output argument)
%clrscl('rmbcgyr',36);
%-------------------
% he tried using values less than 1 so that yellow might be more
% visible on the screen.
q=1/sqrt(2);

cmap=[];
for ic=1:length(pal)
  switch pal(ic);
    case 'r'
      cmap=[cmap; 1 0 0];
    case 'g'
      cmap=[cmap; 0 1 0];
    case 'b'
      cmap=[cmap; 0 0 1];
    case 'C'
      cmap=[cmap; 0 q q];
    case 'c'
      cmap=[cmap; 0 1 1];
    case 'Y'
      cmap=[cmap; q q 0];
    case 'y'
      cmap=[cmap; 1 1 0];
    case 'M'
      cmap=[cmap; q 0 q];
    case 'm'
      cmap=[cmap; 1 0 1];
    case 'w'
      cmap=[cmap; 1 1 1];
    case 'k'
      cmap=[cmap; 0 0 0];
    otherwise
      error('palette error');
  end
end

nmap=size(cmap,1);
if nmap<2
  error('palette error');
end

if nargin<2
  ncol=nmap;
end

ncol=floor(ncol);
if ncol<nmap
  error('palette error');
end

xi=[0:nmap-1]/(nmap-1);
xo=[0:ncol-1]/(ncol-1);

cmap=interp1(xi,cmap,xo,'linear');

if nargout>0
  return
end
% if no arguments, plot the color values, and show the
% palette.
disp('no output args; just plot color map')

xbar=[1:ncol];

clf
plot( xbar, cmap(:,1), 'r');
hold on
plot( xbar, cmap(:,1), 'ro');
plot( xbar, cmap(:,2), 'g');
plot( xbar, cmap(:,2), 'gx');
plot( xbar, cmap(:,3), 'b');
plot( xbar, cmap(:,3), 'b^');
hold off
grid on;
set( gca, 'xlim', [1 ncol]);

colormap(cmap);
colorbar

disp('hit key to continue')
pause

% show the palette again, for an arbitrary full scale
% threshold to 0. (based on his db scales)
threshold=-30; ntv=-threshold/ncol;

scale=[threshold:(0-threshold)/(ncol):0]-ntv/2;
clf
imagesc( scale, [0 1]', [scale' scale']', [threshold 0]);
set( gca, 'position',  [0.1 0.1 0.8 0.8]);
hold on
contour( scale, [0 1]', [scale' scale']', [threshold:3:0], 'k')
hold off
set( gca, 'xtick', [threshold 0]);
set( gca, 'ytick', []);
axis([ threshold 0 0 1]);
xlabel('Power (dB)');

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
