if (exist('inputfil')==0)
   filenm=input('SAGA filename ? ','s')
else
   filenm=inputfil
end;
eval( filenm);
iforward=iopt(30);
eval(['load '  filenm '.ext -ascii']);
eval(['fit='  filenm '; clear ' filenm ])
figure(3)
clf
axs = [0,128, 0.3 .9];
  t1 = fit(:,3);
  bart=1-t1;
  [maxbart,index]= max(bart);

for j=1:nparm
  t2 =fit(:,j+3);
  for jj=0:npdpsamp(j)-1;
    xsample=find(t2<=jj+2 &  t2>=jj-2 );
    for ibart=1:100
      ybart(ibart)=-0.005+ibart*0.01;
      yup=ybart(ibart)+0.005; ylow=ybart(ibart)-0.015;
      temp(ibart,jj+1)=size(find(bart(xsample)<yup & bart(xsample)>ylow),1);
    end
  end
  numbart(:,:,j)=temp;
end

maxamp=max(max(max( numbart)));
cutup=maxamp/5

for j=1:nparm
  subplot(ceil(nparm/2),2,j)
  t2 =fit(:,j+3);
  l1=f_min(j); l2=f_max(j);
  if par2phy(j)==9 
     t2=t2/1000;
  end
  dsamp=(l2-l1)/(npdpsamp(j)-1);
  xval=l1+[0:(npdpsamp(j)-1)]*dsamp;
%   plot(0,0)
%  plot(t2,1-t1,'k.','MarkerSize',1)
%  hold on
val=hist(t2,128);
val=filter([1 1 1]/3, 1 , val);
val=val/max(val);
%plot(val)
for jj=0:127
  jj;
  xsample=find(t2<=xval(jj+1)+2*dsamp &  t2>=xval(jj+1)-2*dsamp );
  if (xsample>0) 
     enve(jj+1)=max(bart(xsample));
  else
     enve(jj+1)=0;
  end
%  for ibart=1:100
%   ybart(ibart)=-0.005+ibart*0.01;
%   yup=ybart(ibart)+0.005; ylow=ybart(ibart)-0.015;
%   find(bart(xsample)<yup & bart(xsample)>ylow)
%   numbart(ibart,jj+1)=size(find(bart(xsample)<yup & bart(xsample)>ylow),1);
%  end
end
imagesc(xval,ybart,numbart(:,:,j),[0 cutup]) %,colorbar
hold on
axis('xy')

plot(xval,enve,'r')   % plot(filter([1 1 1 ]/3,1,val),'g')

  plot([(t2(index)) (t2(index))],[1 min([maxbart+0.03 1])],'Color',[.85,.85,.85],'Linewidth',2); 
axis([l1 l2 axs(3) axs(4)])
% set(gca,'XTick',[0:32:128])
%x = linspace(l1,l2,5);
%set(gca,'XTickLabel',x);
set(gca,'YTick',[axs(3):.1:axs(4)])
set(gca,'Box','off')

xlabel( xtitles(iforward,par2phy(j)),'Fontsize',10);
end

colormap(flipud(gray));


 axes('position',[0,0,1,1]); axis('off')
 text(.05,.02, ['File: ',pwd,'/',filenm,'; Date: ',date],'fontsize',6)
orient('tall')