%field.m test algorithm for 3D array rotations
% variables, tau, eta(1:3),
%
% z axis is down, x axis away from source, y axis to the rigth of x axis
% this is a right handed coordinate system
%tau=maximum displacement between curved and linear array with
%	same endpoints
%
%eta(3)-swivel  (counterclockwise around z-axis )
%eta(1)-tilt    (counterclockwise around Y-axis ) 
%eta(2)-azimuth (counterclockwise around z-axis )

tau=-25;  10;
eta(1)=90;    15;        % rotation 1
eta(2)=90;    45; 45;    % rotation 2
eta(3)=110;    45;  30;   % rotation 3


z=[
  0.00
3.4   
7.14  
11.32 
15.94 
26.26 
32.17 
46.5  
54.95 
64.18 
77.17 
85.03 
96.9  
107.53
151.40
154.82
158.57
167.08
172.15
177.77
184.06
199.48
206.75
215.98
225.99
236.69
259.85
]';
I=1:length(z);Nel=length(z);
x=zeros(size(z));y=x;
rngref=10000;

%First, subtract water depth from z
%z=D-z;
%Actually, subtract top phone depth
zref=z(1);
z=-zref+z;

zref=0
%Generate array curvature-- Aaron
%x=-4*tau*(I-1)./(Nel-1)+4*tau*((I-1).^2/((Nel-1).^2));
% Peter: changed sign; pos tau => pos x;
%        z coordinate itroduced otherwice it will not work for non
%        uniform array
x= 4*tau* ( z./z(Nel)- (z./z(Nel)).^2 ) ;
len=abs(z(1)-z(Nel));
if (len<tau) 
   disp(' This is only a small correction') 
end;
% shorten array using first order:
z=z*(1-(3/8*tau^2/len^2) );

% Swivel around z axis
[xp,yp]    =rot(x,y,eta(1)); zp=z;
% tilt
% wrong [xpp,zpp]  =rot(xp,zp,eta(1))    ;ypp=yp;
[zpp,xpp]  =rot(zp,xp,eta(2))    ;ypp=yp;
% azimuth
[xppp,yppp]=rot(xpp,ypp,eta(3)); zppp=zpp;

r=sqrt((rngref+xppp).^2+yppp.^2);

z   =z   +zref;
zp  =zp  +zref;
zpp =zpp +zref;
zppp=zppp+zref;

figure
subplot(2,2,1)
plot3(x,y,z,'o');
%fill3(x,y,z,'r');
hold on; plot3([x(1) x(Nel)] ,[y(1) y(Nel)],[z(1) z(Nel)],'r');
xlabel('x');ylabel('y');zlabel('z')

a(1,:)=axis;
title(['a) Parabola, bow=' ,num2str((tau))], 'fontsize',12 )

subplot(2,2,2)
plot3(xp,yp,zp,'o');
%fill3(xp,yp,zp,'r');
hold on; plot3([xp(1) xp(Nel)] ,[yp(1) yp(Nel)],[zp(1) zp(Nel)],'r');
xlabel('x');ylabel('y');zlabel('z')
title(['b) Rot1, \theta_1=' ,num2str(eta(1))], 'fontsize',12 )
a(2,:)=axis;

subplot(2,2,3)
plot3(xpp,ypp,zpp,'o');
%fill3(xpp,ypp,zpp,'r');
hold on; plot3([xpp(1) xpp(Nel)] ,[ypp(1) ypp(Nel)],[zpp(1) zpp(Nel)],'r');
xlabel('x');ylabel('y');zlabel('z')
title(['c) Rot2, \theta_2=' ,num2str(eta(2))], 'fontsize',12 )
a(3,:)=axis;

subplot(2,2,4)
hand=plot3(xppp,yppp,zppp,'o','LineWidth',2,...);
                'MarkerEdgeColor','k','MarkerFaceColor',[.49 1 .63])
%                    'MarkerSize',12)
plot3(xppp,yppp,zppp,'o');
%fill3(xppp,yppp,zppp,'r');
hold on; plot3([xppp(1) xppp(Nel)] ,[yppp(1) yppp(Nel)],[zppp(1) zppp(Nel)],'r');
xlabel('x');ylabel('y');zlabel('z')
title(['d) Rot3, \theta_3=' ,num2str(eta(3))], 'fontsize',12 )
a(4,:)=axis;

ax1=max(a); ax2=min(a);
ax=[ax2(1) ax1(2) ax2(3) ax1(4) ax2(5) ax1(6)];
ax=[-80 256 -50 256  0 256];

subplot(2,2,1)
axis(ax); set(gca,'zdir','rev','ydir','rev')
grid, box

subplot(2,2,2)
axis(ax); set(gca,'zdir','rev','ydir','rev')
grid, box
plot3( [0 0],[0 0], [ax(5) ax(6)],'go-','linewidth',1);

subplot(2,2,3)
axis(ax); set(gca,'zdir','rev','ydir','rev')
grid, box
plot3( [0 0],[ax(3) ax(4)], [zref  zref],'go-','linewidth',1);

subplot(2,2,4)
axis(ax); set(gca,'zdir','rev','ydir','rev')
%axis('equal');grid
grid, box
plot3( [0 0],[0 0], [ax(5) ax(6)],'go-','linewidth',1);

% directional angles
%di=[(xpp(Nel)-xpp(1)) (ypp(Nel)-ypp(1)) (zpp(Nel)-zpp(1))];
di=[(xppp(Nel)-xppp(1)) (yppp(Nel)-yppp(1)) (zppp(Nel)-zppp(1))];
dist=(di*di')^.5;
di=di/dist;

angle=acos(di)*180/pi;

di=[(xppp-xppp(1))' (yppp-yppp(1))' (zppp-zppp(1))'];



%for jj=2:Nel
% dist2=(di(jj,:)*di(jj,:)')^.5;
% cr=cross(di(jj,:),di(Nel,:))/dist2/dist;
% cr=cr/(cr*cr')^.5;
% angle=acos(cr)*180/pi;
%end
