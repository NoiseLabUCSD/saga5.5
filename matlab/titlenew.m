function h1 = titlenew(string,varargin);
%function h1 = titlenew(string,horizonalig,titlefont);
% to justify the title text to appear on the top left-hand or
% right-hand side of the figure
% default is 'left';
% cfh, 1/5/2005   varargin{n}(2)
   titlefont = 14;
   horizonalig = 'left';

if nargin < 1
   error('Requires at last 1 input argument.')
elseif (nargin == 2)
   horizonalig = 'left';
   titlefont = 14;
   if isstr(varargin{:}) 
      horizonalig = varargin{:};
   end
   if ~isstr(varargin{:})
      titlefont = varargin{:};
   end
elseif nargin == 3
   for tt = 1:2
      if isstr(varargin{tt}) 
         horizonalig = varargin{tt};
      end
      if ~isstr(varargin{tt})
         titlefont = varargin{tt};
      end
   end
end
t = title(string,'fontsize',titlefont);

set(t, 'horizontalAlignment', horizonalig);
set(t, 'units', 'normalized');
h1 = get(t, 'position');
   
%   keyboard
if strcmp(horizonalig,'left');
   set(t, 'position', [0.01 h1(2)*0.99]);
elseif strcmp(horizonalig,'right');
   set(t, 'position', [1-0.01 h1(2) h1(3)]);
end

