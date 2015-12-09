function []=wplot(arg1,arg2,arg3,arg4);
% WPLOT        Creates wiggle trace, with
%	  areas above zero filled.
%	  (works like plot, exept that
%          only one function can be potted
%          at the time)
%     
%	  syntax: wplot(x,'c-');
%             or:      (x,y,'c-');
%             or:      (x,y,offs,'c-');
%         offs  (default=0)
%         'c-'  color and linetype (optional)  
%	  5/1/95 Taco van der Leij
%    Peter Gerstoft 8/97: The number of samples has to bee dense; otherwise
%            the filling does not work around zero-crossings. 
%        [Fix:  The x-coordinates must be recalculated when there is a
%            zero crossing.]
nargin_wplot=nargin;
val=eval(['arg' num2str(nargin_wplot)]);
if isstr(val),
	nargin_wplot=nargin_wplot-1;
	arg=val;
else
	arg='y-';
end;
if nargin_wplot<3,
	offs=0;
else
        offs=arg3;
end;
if nargin_wplot==1,
	y_ax=arg1;
	x_ax=1:length(y_ax);
else
	x_ax=arg1;
	y_ax=arg2;
end;

np=get(gca,'nextplot');
nx=length(x_ax);
ny=length(y_ax);
if size(x_ax,1)==nx,
	x_ax=x_ax.';
end;
if size(y_ax,1)==ny,
	y_ax=y_ax.';
end;
sx=sort(x_ax);
plot(x_ax,y_ax,arg);
%keyboard
set(gca,'nextplot','add');

if (sx==x_ax)|(sx==x_ax(nx:-1:1)), 
	ym=offs;
	dy=y_ax-ym;
%keyboard
	patch([x_ax(1) x_ax x_ax(nx)],...
	     [ym ym+(abs(dy)+dy)/2 ym],arg);
else,
	xm=offs;
	dx=x_ax-xm;
	fill([xm xm+(abs(dx)+dx)/2 xm],...
	     [y_ax(1) y_ax y_ax(ny)],arg);
end;
set(gca,'nextplot','add');
plot(x_ax,y_ax,arg);
set(gca,'nextplot',np);





