filenm='asiaexcsdm135.trf'

%[s,par]        = readascii(filenm)
[s1,par]        = rereadtrf(filenm);
  ranges=par(5)+par(6)*( 0:par(7)-1);
  depth= par(2)+ (par(3)-par(2))/(par(4)-1)*(0:par(4)-1);
s1=-s1;  % fit to TL !
  
ntrf=size(s1,3);
nf=1; ifreq=1;
nr=par(7)
nd=par(4)
%for ii=1:ntrf
%  tlcov{:,:ii}=s1{ii};
%end

figure
imagesc( ranges,depth,squeeze(s1(:,:,1))),colorbar


figure
[xval,marg,fitval]=  plotenum('asiaexcsdm135',3);
xlabel('Bottom sound speed (m/s)')
expfit = exp(-fitval);
expfit=expfit(1:ntrf);
lss=expfit/sum(expfit);
%print -depsc marg2d



xval=40:65
id=nd/2
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

figure
subplot(2,1,1)
imagesc(ranges,xval,pa',cbound),colorbar,axis('xy')
xlabel('Range (m)'), ylabel('Transmission loss (dB)'); style1
title('Prior h = 50 m')
hold on; axis ij
plot([830 830],[40 65],'w','linewidth',2)
subplot(2,1,2)
imagesc(ranges,xval,pp',cbound),colorbar,axis('xy')
xlabel('Range (m)'), ylabel('Transmission loss (dB)'); style1
title('Posterior h = 50 m')
hold on; axis ij
plot([830 830],[40 65],'w','linewidth',2)
%print -depsc  postrange


figure
subplot(2,1,1)
plot(xval,pp(83,:),'b',xval,pa(80,:),'b--','linewidth',1)
ylabel('Probability')
xlabel('Transmission loss (dB)'); style1
print -depsc tl50

figure
cbound=[30 90]
plot(squeeze(s1(id,:,:)))
set(gca,'ydir','rev')

for id=1:nd
for ir=1:nr
     ss=squeeze(s1(id,ir,:));
     [no,xo]= mlpdf(ss,lss,xval);
     ss1=0; ii1=0; ii2=0; ii3=0;
     for ii=1:length(xval)
        ss1=ss1+no(ii);
        if (ss1>0.05 &ii1==0); ii1=ii; end
        if (ss1>0.95 & ii2==0); ii2=ii; end	
        if (ss1>0.5 & ii3==0); ii3=ii; end	
     end
     spread(ir,id)=xval(ii2)-xval(ii1);    med(ir,id)=xval(ii3);
     lval(ir,id)=xval(ii1); uval(ir,id)=xval(ii2);
     %
     [no1,xo]= mpdf(ss,xval);
     ss1=0; ii1=0; ii2=0; ii3=0;
     for ii=1:length(xval)
        ss1=ss1+no1(ii);
        if (ss1>0.05 &ii1==0); ii1=ii; end
        if (ss1>0.95 & ii2==0); ii2=ii; end	
        if (ss1>0.5 & ii3==0); ii3=ii; end	
     end
     spreadp(ir,id)=xval(ii2)-xval(ii1);   medp(ir,id)=xval(ii3);
     lvalp(ir,id)=xval(ii1); uvalp(ir,id)=xval(ii2);
end
end


figure
cbound=[35 65]
subplot(2,1,1)
imagesc(ranges,depth,medp',cbound),colorbar

xlabel('Range (m)'), ylabel('Depth (m)'); style1
title('Prior median')

subplot(2,1,2)
imagesc(ranges,depth,med',cbound),colorbar
xlabel('Range (m)'), ylabel('Depth (m)'); style1
title('Posterior median')
colormap(fliplr(jet))
print -depsc mediantl
figure
cbound=[0 20]
subplot(2,1,1)
imagesc(ranges,depth,spreadp',cbound),colorbar

%xlabel('Range (m)'), ylabel('Depth (m)'); style1
title('a) Prior, range of 5th to 95th percentiles')
style1
subplot(2,1,2)
imagesc(ranges,depth,spread',cbound),colorbar

xlabel('Range (m)'), ylabel('Depth (m)'); style1
title('b) Posterior, range of 20th to 80th percentiles')
style1
colormap(jet)


print -depsc spread

figure
cbound=[0 25]
subplot(2,1,1)
imagesc(ranges,depth,spreadp',cbound),colorbar

%xlabel('Range (m)'), ylabel('Depth (m)'); style1
title('a) Prior, range of 5th to 95th percentiles')

subplot(2,1,2)
imagesc(ranges,depth,spread',cbound),colorbar

xlabel('Range (m)'), ylabel('Depth (m)'); style1
title('b) Posterior, range of 5th to 95th percentiles')
colormap(jet)


print -depsc spread

id=nd/2
figure
subplot(2,1,1)
fill([ranges fliplr(ranges)],[uvalp(:,id)' fliplr(lvalp(:,id)')],[0.9 0.9 ...
      0.9])
hold on;
plot(ranges,med(:,id),'r-') %,ranges,tc1tlref(id,:)','k-')
set(gca,'ydir','rev','ylim',[40 80])
title(['prior, depth=',num2str(depth(id))])
style1, text(10,37,'a)','fontsize',14)
%
subplot(2,1,2)
fill([ranges fliplr(ranges)],[uval(:,id)' fliplr(lval(:,id)')],[0.9 0.9 ...
      0.9])
hold on;
plot(ranges,med(:,id),'r-') %,ranges,tc1tlref(id,:)','k-')
%plot(ranges,uval(:,id),'b--',ranges,lval(:,id),'b--',ranges,med(:,id),'k-')
set(gca,'ydir','rev','ylim',[40 80])
title(['posterior, depth=',num2str(depth(id))])
  xlabel('Ranges (m)'); ylabel('Tloss (dB)')
style1, text(10,37,'b)','fontsize',14)
% 
%mean(abs(med(:,id)-tc1tlref(id,:)'))
%mean(abs(medp(:,id)-tc1tlref(id,:)'))
print -depsc asiaextlvar50
