function [no,xo] = mlpdf(y,lh,x)
%HIST  Histogram.
%   N = HIST(Y,lh) bins the elements of Y into 10 equally spaced containers
%   and returns the number of elements in each container.  If Y is a
%   matrix, HIST works down the columns.
%
%   N = HIST(Y,lh,M), where M is a scalar, uses M bins.
%
%   N = HIST(Y,X), where X is a vector, returns the distribution of Y
%   among bins with centers specified by X.  Note: Use HISTC if it is
%   more natural to specify bin edges instead.
%
%   [N,X] = HIST(...) also returns the position of the bin centers in X.
%
%   HIST(...) without output arguments produces a histogram bar plot of
%   the results.
%
%   See also HISTC.

%   J.N. Little 2-06-86
%   Revised 10-29-87, 12-29-88 LS
%   Revised 8-13-91 by cmt, 2-3-92 by ls.
%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 5.20 $  $Date: 2002/06/05 17:06:38 $

if nargin <1
    error('Requires 2 or 3 input arguments.')
end
if nargin == 2
    x = 10;
end
[m,n] = size(y);
if (size(lh,1)~=m | size(lh,2)~=n)
  error('x and likelihood must have identical size')
end
if min(size(y))==1, lh=lh(:); y = y(:); end
if isstr(x) | isstr(y)
    error('Input arguments must be numeric.')
end

[m,n] = size(y);
if isempty(y),
    if length(x) == 1,
       x = 1:x;
    end
    nn = zeros(size(x)); % No elements to count
else
    if length(x) == 1
        miny = min(min(y));
        maxy = max(max(y));
    	  if miny == maxy,
    		  miny = miny - floor(x/2) - 0.5; 
    		  maxy = maxy + ceil(x/2) - 0.5;
     	  end
        binwidth = (maxy - miny) ./ x;
        xx = miny + binwidth*(0:x);
        xx(length(xx)) = maxy;
        x = xx(1:length(xx)-1) + binwidth/2;
    else
        xx = x(:)';
        miny = min(min(y));
        maxy = max(max(y));
        binwidth = [diff(xx) 0];
        xx = [xx(1)-binwidth(1)/2 xx+binwidth/2];
        xx(1) = min(xx(1),miny);
        xx(end) = max(xx(end),maxy);
    end
    nbin = length(xx);
    % Shift bins so the internal is ( ] instead of [ ).
    xx = full(real(xx)); y = full(real(y)); % For compatibility
    bins = xx + max(eps,eps*abs(xx));
%%%%%    nn = histc(y,[-inf bins],1);
%%%%%    keyboard
    pgbin=[-inf bins inf];
    for I=1:length(pgbin)-1
      [ii,aa]= find(y    <pgbin(I+1));
      [jj,aa]= find( y(ii)>pgbin(I));
      nn(I,:)=sum(lh(ii(jj)));
    end
%    keyboard
  % Combine first bin with 2nd bin and last bin with next to last bin
    nn(2,:) = nn(2,:)+nn(1,:);
    nn(end-1,:) = nn(end-1,:)+nn(end,:);
    nn = nn(2:end-1,:);
end
    nn=nn/sum(nn); %pg norm to probablility

if nargout == 0
    bar(x,nn,'hist');
else
  if min(size(y))==1, % Return row vectors if possible.
    no = nn';
    xo = x;
  else
    no = nn;
    xo = x';
  end
end
