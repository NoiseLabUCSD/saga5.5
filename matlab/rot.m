function [x1p,x2p]=rot(x1,x2,thta)
	thta=thta*pi/180;
	R=[cos(thta) -sin(thta);
	   sin(thta) cos(thta)];
	out=R*[x1;x2];
	x1p=out(1,:);
	x2p=out(2,:);
	