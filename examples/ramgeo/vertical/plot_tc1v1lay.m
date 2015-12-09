inputfil='tc1_v1j'
%inputfil='tc1_v6j'
inputfil='tc1_v6j'
inputfil='tc1v6lay'
inputfil='tc1v1lay'

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
xscale=8000;    xoff=300;
for j1=0:4:0   %4:16
  x1=xscale*(a(:,j1+1)-mean(a(:,j1+1)))+xoff;
  x2=xscale*(a(:,j1+3)-mean(a(:,j1+3)))+xoff;
plot(x1,a(:,j1+2),'k-',x2,a(:,j1+4),'r','linewidth',2)
hold on
  plot([xoff xoff],[20 80],'b--')
  xoff=xoff+250;
end
set(gca,'Ydir','rev','Ylim',[20 80],'Xtick',[50 300])
subplot(4,6,19:20); hold off

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
xx={'sound speed s1 (m/s)', 'sound speed incr s2 (m/s)', 'sound speed incr s3 (m/s)', 'sound speed incr @bot+ (m/s)', ...
      'Thickness (m)','Thickness (m)','Thickness (m)','Density (g/cm3)','Density (g/cm3)', 'Attenuation (dB/\lambda)', 'Attenuation (dB/\lambda)','Bathymetry (m)', 'Bathymetry (m)'};

text(-0.2 ,   0.35,xx   ,'Fontsize',7);
  text(0.6  , 0.35,xbest,'Fontsize',7);
  text(0.8  , 0.35,xstd ,'Fontsize',7);
axis off


subplot(4,6,21:22); hold off
ssp=res(1,1);   %0
ssp(2)=ssp(1)
ssp(3)=ssp(1)+res(2,1);  %3
ssp(4)=ssp(3)
ssp(5)=ssp(3)+res(3,1); %10
ssp(6)=ssp(5)
ssp(7)=ssp(6)+res(4,1); %35
zdep=0
zdep(2)=zdep(1)+res(5,1) %3
zdep(3)=zdep(2)
zdep(4)=zdep(3)+res(6,1) %10
zdep(5)=zdep(4)
zdep(6)=zdep(5)+res(7,1) %35
zdep(7)=zdep(6)

plot([ssp ssp(7)],[zdep zdep(7)+5],'r')
hold on



%plot([1525 1610 1800 1800],[0   35 35 res(6)+5],'k')
xlabel(' ssp (m/s)')
ylabel(' depth (m)')
axis([1500 max([ssp(7)+40 1900]) 0  max([zdep(7)+5 40])]); set(gca,'Ydir','rev')


subplot(2,3,6)
axis off

subplot(2,3,6)
axis off

tmpjpg=imread('temp.jpg','jpeg');
imagesc(tmpjpg)
axis off
title(['best fit:  ' num2str(bestfit) '  '])

axes('position',[0,0,1,1]); axis('off')
text(.05,.01, ['File: ',pwd,'/',inputfil,'; Date: ',date],'fontsize',6)

eval(['print -depsc ' filenm]) 
eval(['print -dpng ' filenm]) 





