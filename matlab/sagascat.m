if (exist('inputfil')==0)
   filenm=input('SAGA filename ? ','s')
else
   filenm=inputfil;
end;
eval( filenm);
iforward=iopt(30);
eval(['load '  filenm '.ext -ascii']);
eval(['fit='  filenm '; clear ' filenm ]);
figure
clf
axs = [0,128, 0.3 .9];
  t1 = fit(:,3);
  bart=1-t1;
  [maxbart,index]= max(bart);

  nparm=size(fit,2)-3

t1db=-10*log10(t1);
t1min=min(t1db);t1max=max([0 max(t1db)+0.2]);

if (length(t1)<100);
 msize=8;
elseif (length(t1)<1000);
 msize=6;
elseif (length(t1)<5000);
 msize=3;
elseif (length(t1)<20000);
 msize=2;
end

for j=1:nparm
  subplot(ceil(nparm/2),2,j)
  t2 =fit(:,j+3);
  l1=f_min(j); l2=f_max(j);
  if par2phy(j)==9; 
     t2=t2/1000;
  end
  dsamp=(l2-l1)/(npdpsamp(j)-1);
  xval=l1+[0:(npdpsamp(j)-1)]*dsamp;
%   plot(0,0)
%  plot(t2,1-t1,'k.','MarkerSize',1)
%  hold on
plot(t2,t1db,'k.','markersize',6);
axis([l1 l2 t1min-1 t1max]);
xlabel( xtitles(iforward,par2phy(j)),'Fontsize',10);
ylabel('- 10 log(\phi)')
end

%nhq 20140902
% axes('position',[0,0,1,1]); axis('off')
% text(.05,.01, ['File: ',pwd,'/',filenm,'; Date: ',date],'fontsize',6)
orient('tall')
