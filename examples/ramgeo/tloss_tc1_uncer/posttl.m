% the matlab file backup 

filenm='tc1v1lay'

%rereadramgeo
mm=10000
filenm='tc1v1layresa'
   [s1,par]  = rereadtrf([filenm '.trf'],mm);
   ranges = par(5)+par(6)*(0:par(7)-1);
   depth  = par(2)+(par(3)-par(2))/(par(4)-1)*(0:par(4)-1);
   s1 = -s1;  % fit to TL !
%[s,par]        = readascii(filenm)
%[s1,par]        = rereadtrf(filenm);
 % ranges=25*(0:nx-1)
 % depth= 2*2*(0:nz/2)
%s1=-s1;  % fit to TL !
% s1=s;
 nx=par(7)
 nr=nx;
 md=par(4)
 nd=par(4)
ntrf=size(s1,3);
nf=1; ifreq=1;
%nr=nx
%nd=nz/2+0.5
%for ii=1:ntrf
%  tlcov{:,:ii}=s1{ii};
%end

%figure
%imagesc( ranges,depth,squeeze(s1(:,:,1))),colorbar


figure
filenm='tc1v1lay'
[xval,marg,fitval]=  plotenum(filenm,3);
subplot(4,4,1), xlabel('Layer-1 sound speed (m/s)'),set(gca,'YTickLabel',[])
subplot(4,4,6), xlabel('Layer-2 sound speed incr (m/s)'),set(gca,'YTickLabel',[])
subplot(4,4,11), xlabel('Sound speed incr (m/s)'),set(gca,'YTickLabel',[])
subplot(4,4,16), xlabel('Sound speed incr (m/s)'),set(gca,'YTickLabel',[])
expfit = exp(-fitval);
expfit=expfit(1:ntrf);
lss=expfit/sum(expfit);
colormap('gray');colormap(1-colormap)

%print -depsc tc1marg2d
%%

xval=40:80;
id=1; %nd/2
for ir=1:nr
     ss=squeeze(s1(id,ir,:));
     [no,xo]= mlpdf(ss,lss,xval);
     pp(ir,:)=no;
     [no1,xo]= mpdf(ss,xval);
     pa(ir,:)=no1;       
 end

%end
cbound=[0 0.15];
%imagesc(ranges,height,pp',cbound),colorbar,axis('xy')
ir=82;
figure
subplot(2,1,1)
imagesc(ranges,xval,pa',cbound),colorbar,axis('xy')
xlabel('Range (m)'), ylabel('Transmission loss (dB)'); style1
%title('Prior')
text(100,36,'(a)','fontsize',12)
hold on; axis ij
plot([ranges(ir) ranges(ir)],[40 80],'k--','linewidth',2)
%%
subplot(2,1,2)
imagesc(ranges,xval,pp',cbound),colorbar,axis('xy')
xlabel('Range (m)'), ylabel('Transmission loss (dB)'); style1
%title('Posterior')
hold on; axis ij
plot([ranges(ir) ranges(ir)],[40 80],'k--','linewidth',2)
text(100,36,'(b)','fontsize',12)
%plot([830 830],[40 65],'w','linewidth',2)
%print -deps  tc1postrange
%colormap('gray');colormap(1-colormap)
%4
%%
ir=80;
figure
subplot(2,1,1)
plot(xval,pp(ir,:),'b',xval,pa(ir,:),'b--','linewidth',1)
ylabel('Probability')
xlabel('Transmission loss (dB)'); style1
axis ([40 80 0 0.6])
print -depsc tl42
%%
figure
cbound=[30 90]
plot(squeeze(s1(id,:,:)))
set(gca,'ydir','rev')

for id=1:nd
for ir=1:nr
     ss=squeeze(s1(id,ir,:));
     [no,xo]= mlpdf(ss,lss,xval);
     sp=spread(no,xo);
     spreadarr(ir,id)=sp(2)-sp(1);    med(ir,id)=sp(3);
     lval(ir,id)=sp(1); uval(ir,id)=sp(2);
     %
     [no1,xo]= mpdf(ss,xval);
     sp=spread(no1,xo);
     spreadarrp(ir,id)=sp(2)-sp(1);   medp(ir,id)=sp(3);
     lvalp(ir,id)=sp(1); uvalp(ir,id)=sp(2);
end
end

figure
cbound=[35 65]
subplot(2,1,1)
imagesc(ranges,depth,medp',cbound),colorbar
hold on; plot([0 5000],[90 150],'k','linewidth',2)
xlabel('Range (m)'), ylabel('Depth (m)'); style1
colormap('gray')
%title('Prior median')
text(100,-10,'(a)','fontsize',12)

subplot(2,1,2)
imagesc(ranges,depth,med',cbound),colorbar
xlabel('Range (m)'), ylabel('Depth (m)'); style1
hold on; plot([0 5000],[90 150],'k','linewidth',2)
%title('Posterior median')
text(100,-10,'(b)','fontsize',12)
colormap(fliplr(jet))
colormap('gray')
print -deps mediantl


figure
cbound=[5 20]
subplot(2,1,1)
imagesc(ranges,depth,spreadarrp',cbound),colorbar
colormap('gray');colormap(1-colormap)
hold on; plot([0 5000],[90 150],'k','linewidth',2)
xlabel('Range (m)'), ylabel('Depth (m)'); style1
text(100,-10,'(a)','fontsize',12)

%xlabel('Range (m)'), ylabel('Depth (m)'); style1
%title('a) Prior, range of 5th to 95th percentiles')
style1
subplot(2,1,2)
imagesc(ranges,depth,spreadarr',cbound),colorbar
hold on; plot([0 5000],[90 150],'k','linewidth',2)

%title('b) Posterior, range of 5th to 95th percentiles')
xlabel('Range (m)'), ylabel('Depth (m)'); style1
text(100,-10,'(b)','fontsize',12)
%colormap(jet)
colormap('gray');colormap(1-colormap)

print -deps spread


Keyboard
STOP
FILE NEED SOME TESTING BELOW IN ORDER TO WORK
%%%%%%%%%%%%%%%%%%%%%%%%%%%single diagramid=21; ...
%nd/2
%%
figure
for id=1:nd
subplot(2,1,1)
plot(ranges,uvalp(:,id),'b--',ranges,lvalp(:,id),'b--',ranges,medp(:,id),'k-')
set(gca,'ydir','rev')
title(['prior, depth=',num2str(depth(id))])
subplot(2,1,2)
plot(ranges,uval(:,id),'b--',ranges,lval(:,id),'b--',ranges,med(:,id),'k-')
set(gca,'ydir','rev')
title(['prior, depth=',num2str(depth(id))])
  xlabel('Ranges (m)'); ylabel('Tloss (dB0)')
  pause
end

load tc1tlref
%  hold on; plot(ranges,tc1tlref(id,:)','r:')
id=27
subplot(2,1,1)
fill([ranges fliplr(ranges)],[uvalp(:,id)' fliplr(lvalp(:,id)')],[0.9 0.9 ...
      0.9])
hold on;
plot(ranges,medp(:,id),'r-',ranges,tc1tlref(id,:)','k-')
set(gca,'ydir','rev','ylim',[40 80])
title(['prior, depth=',num2str(depth(id))])
style1, text(10,37,'a)','fontsize',14)
%
subplot(2,1,2)
fill([ranges fliplr(ranges)],[uval(:,id)' fliplr(lval(:,id)')],[0.9 0.9 ...
      0.9])
hold on;
plot(ranges,med(:,id),'r-',ranges,tc1tlref(id,:)','k-')
%plot(ranges,uval(:,id),'b--',ranges,lval(:,id),'b--',ranges,med(:,id),'k-')
set(gca,'ydir','rev','ylim',[40 80])
title(['posterior, depth=',num2str(depth(id))])
  xlabel('Ranges (m)'); ylabel('Tloss (dB)')
style1, text(10,37,'b)','fontsize',14)
 
print -depsc tc1tlvar104
%%
id=21
subplot(2,1,1)
fill([ranges fliplr(ranges)],[uvalp(:,id)' fliplr(lvalp(:,id)')],[0.9 0.9 ...
      0.9])
hold on;
plot(ranges,medp(:,id),'k--',ranges,tc1tlref(id,:)','k-')
set(gca,'ydir','rev','ylim',[40 80])
%title(['prior, depth=',num2str(depth(id))])
 xlabel('Range (m)'); ylabel('Transmission loss (dB)');
style1, 
text(10,37,'(a)','fontsize',14)
%
subplot(2,1,2)
fill([ranges fliplr(ranges)],[uval(:,id)' fliplr(lval(:,id)')],[0.9 0.9 ...
      0.9])
hold on;
plot(ranges,med(:,id),'k--',ranges,tc1tlref(id,:)','k-')
%plot(ranges,uval(:,id),'b--',ranges,lval(:,id),'b--',ranges,med(:,id),'k-')
set(gca,'ydir','rev','ylim',[40 80])
%title(['posterior, depth=',num2str(depth(id))])
  xlabel('Range (m)'); ylabel('Transmission loss (dB)');
style1, 
text(10,37,'(b)','fontsize',14)
% 
sqrt(mean(( med(:,id)-tc1tlref(id,:)').^2))
sqrt(mean((medp(:,id)-tc1tlref(id,:)').^2))
print -deps tc1tlvar80

fill([ranges fliplr(ranges)],[uvalp(:,id)' fliplr(lvalp(:,id)')],[0.9 0.9 0.9])
