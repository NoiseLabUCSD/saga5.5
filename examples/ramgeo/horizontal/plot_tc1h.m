inputfil='tc1_v1j'
inputfil='tc1v1'
%inputfil='tc1_v6j'
inputfil='tc1h2'

plotsaga
print -djpeg temp
%figure(iframe)
iframe=figure

 tc='1'
  casedir=[ '/heart5/gerstoft/itworkshop/tc' tc '/']
ifreq=0
freq=250;
if (freq<100)
  [pr par]  =read_it([casedir 'ivwkt' tc '_v_00' num2str(freq) '.cpr']);
else
  [pr par]    =read_it([casedir 'ivwkt' tc '_v_0' num2str(freq) '.cpr']);
  [prh parh]  =read_it([casedir 'ivwkt' tc '_h_0' num2str(freq) '.cpr']);
end  


depthv=par(5)+par(6)*[0:(par(7)-1)];
rangeh=parh(2)+parh(3)*[0:(parh(4)-1)];
fid=fopen([inputfil '.trf'],'r');
  dr=fscanf(fid,'%f',1); dz=fscanf(fid,'%f',1)
  nr=fscanf(fid,'%d',1);  nz=fscanf(fid,'%d /n',1)
  x=fscanf(fid,'%f',[2*nz nr]);
fclose(fid)
z=x(1:2:2*nz,:)+i*x(2:2:2*nz,:);

range=2*dr:dr:(nr+1)*dr;
depth=dz:dz:nz*dz;

tl=20*log10(abs(z))+ones(size(z,1),1)*10*log10(range/1000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(4,1,1),hold off
imagesc(range,depth,tl,[-75 -45])%,colorbar  
title(['testcase ' inputfil])
hold on
plot([0 5000],[90 150],'k','linewidt',1)
axis([60 5000 0 200])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(4,1,2),hold off
plot(rangeh,20*log10(abs(prh(1,:)))+10*log10(rangeh/1000),'k','linewidth',2)
hold ...
    on
plot(range,tl(25,:),'r','linewidth',1)
plot(rangeh,20*log10(abs(prh(2,:)))+10*log10(rangeh/1000)-25,'k','linewidth',2)
hold on
plot(range,tl(85,:)-25,'r','linewidth',1)
axis([0 5000 -100 -40])

%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%subplot(4,1,3)
subplot(4,3,7:8); hold off
xscale=8000;    xoff=50;
xscale=1;    xoff=0;
for j1=0:4:0   %4:16
  x1=xscale*(a(:,j1+1)-mean(a(:,j1+1)))+xoff;
  x2=xscale*(a(:,j1+3)-mean(a(:,j1+3)))+xoff;
  x1=xscale*(a(:,j1+1))+xoff;
  x2=xscale*(a(:,j1+3))+xoff;
plot(x1,a(:,j1+2)-10*log(x1),'k-','linewidth',2)
hold on
plot(x2,a(:,j1+4)-10*log(x1),'r','linewidth',1)
hold on
%  plot([xoff xoff],[20 80],'b--')
  xoff=xoff+0;   250;
end
set(gca,'Ydir','rev','Ylim',[30 80],'Xtick',[0 3])
set(gca,'Ydir','rev','Ylim',[40 70],'Xtick',[0 3])
subplot(4,6,19:21); hold off

axis off
clear xx,clear xstd, clear xbest
if (res(1,1)>100)
xbest={ int2str(res(1,1))};
else
xbest={ num2str(res(1,1),3)};  
  end
 xstd= {num2str(res(1,2),2)};

 xx = {xtitles(iforward,par2phy(1)) };

for j=2:nparm
%  xt=[xtitles(iforward,par2phy(j)) num2str(res(j,1),4) '  ' ...
%	num2str(res(j,2),3)];
if (res(j,1)>1000)
xbest={char(xbest) int2str(res(j,1))};
else
xbest={char(xbest) num2str(res(j,1),4)};  
  end
xstd= {char(xstd)  num2str(res(j,2),2)};
xx = {char(xx),xtitles(iforward,par2phy(j)) };
end  
xx={'sound speed @0m (m/s)', 'sound speed incr @3m (m/s)', 'sound speed incr @10m (m/s)', 'sound speed incr @bot- (m/s)', 'sound speed incr @bot+ (m/s)', ...
      'sediment thickness (m)','Density (g/cm^3)', 'Attenuation (dB/\lambda)','Bathymetry (m)', 'Bathymetry (m)'};

text(-0.1 , 0.45,xx   ,'Fontsize',8);
  text(0.6  , 0.45,xbest,'Fontsize',8);
  text(0.8  , 0.45,xstd ,'Fontsize',8);
axis off


subplot(4,6,22); hold off
ssp=res(1,1);
ssp(2)=ssp(1)+res(2,1);
ssp(3)=ssp(2)+res(3,1);
ssp(4)=ssp(3)+res(4,1);
ssp(5)=ssp(4)+res(5,1);

plot([ssp ssp(5)],[0 3  10  res(6) res(6) res(6)+5],'r')
hold on
%plot([1525 1610 1800 1800],[0   35 35 res(6)+5],'k')
xlabel(' ssp (m/s)')
ylabel(' depth (m)')
axis([1500 max([ssp(5)+40 1900]) 0  max([res(6)+5 40])]); set(gca,'Ydir','rev')



subplot(2,3,6)
axis off

subplot(2,3,6)
axis off

tmpjpg=imread('temp.jpg','jpeg');
imagesc(tmpjpg)
axis off
title(['best fit:  ' num2str(bestfit) '  '])

axes('position',[0,0,1,1]); axis('off')
text(.05,.02, ['File: ',pwd,'/',inputfil,'; Date: ',date],'fontsize',6)

eval(['print -depsc ' filenm]) 
eval(['print -dpng ' filenm]) 





