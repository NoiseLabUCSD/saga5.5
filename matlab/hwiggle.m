function []=hwiggle(M1,arg1,arg2,arg3,arg4,arg5);
% HWIGGLE       Makes a horizontal wiggle plot of a matrix
%
%         Syntax:
%         HWIGGLE(M,int,offs,scale,fill);
%
%         M    = Seismic record to be plotted
%
%	  If you want to use color or a different
%	  line style, just put the string
%         (e.g. 'r--'), as your last argument
%
%         int   = [ dt dz ]
%         offs  = [ t0 z0 ]
%         scale = scaling factor of the wiggles
%                 use 1 for no overlap
%                (2: default)
%                 if you only use scale as an
%                 argument put it right after M
%         fill  = 0: positive part not filled
%	 	 (1: default)
%
%	  So also: WPLOT CLIP
%         author: Taco van der Leij 9/1/95 

np=get(gca,'nextplot');
nargin_hplot=nargin;
%coline='';
coline='y-';
if nargin_hplot>1,
	val=eval(['arg' num2str(nargin_hplot-1)]);
	if isstr(val),
		nargin_hplot=nargin_hplot-1;
		coline=val;
	end;
end;
if nargin_hplot>1,
	if length(arg1)==1,
		scale=arg1;
		int=[1 1];
	else
		int=arg1;
		if nargin_hplot<4,
			scale=1.5;
		else
			scale=arg3;
		end;
	end;
else
	int=[1 1];
	scale=2;
end;
if (nargin_hplot<3)
	offs=[1 1];
else
	offs=arg2;
end;
if nargin_hplot<5,
	fp=1;
else
	fp=arg4;
end;

dt=int(1);
dx=int(2);
t0=offs(1);
x0=offs(2);
dx2=dx/2;
[nx,nt]=size(M1);

X_ax=(0:dx:(nx-1)*dx)+x0;
T_ax=(0:dt:(nt-1)*dt)+t0;


ud=get(gca,'userdata');
if strcmp(np,'replace')|isempty(ud),
	cla;
	out=max(max(abs(M1)));
	ud=[scale out];
	set(gca,'userdata',ud);
else
	scale=ud(1);
	out=ud(2);
end;

fact=dx2/out*scale;
axis([ t0 t0+(nt-1)*dt ...
       x0-dx2*scale x0+(nx-1)*dx+scale*dx2 ]);

hold on;
%set(gca,'ydir','reverse');
set(gca,'box','on');
if fp==0,
	for x=1:nx, 
		plot(T_ax,X_ax(x)+fact.*M1(x,:),coline)
	end;
else,
	for x=1:nx,
		ddx=fact.*M1(x,:);
		fill([T_ax(1) T_ax T_ax(nt)],...
		     [X_ax(x) X_ax(x)+(abs(ddx)+ddx)/2 X_ax(x)],coline);
		plot(T_ax,X_ax(x)+ddx,'b');
% pg		plot(T_ax,X_ax(x)+ddx,coline);
	end;
end;

set(gca,'nextplot',np);




