clear all
mm = 1000; % mm: # of samples in one fort.81
filen = 'asiaexcsdm135';
nfile = writefort81(filen, mm)
%  which writefort81
%  /mpl/chenfen/saga/matlab/writefort81.m
% create fort.81_...

filenm = [filen 'resa'];
fre = 295;

nfile = 1;
iplot = false; 
iplot = true;

for ij = 1: nfile 
   unix(['/bin/cp -f fort.81_' num2str(ij) ' fort.81']);
   unix(['sed -e s/Fretemp/',num2str(fre),...
	 '/g ',filenm,'.temp > ',filenm,'.dat']);
   
   unix(['resa ' filenm ' snap > out']);
   [s1,par]  = rereadtrf([filenm '.trf'],mm);
   ranges = par(5)+par(6)*(0:par(7)-1);
   depth  = par(2)+(par(3)-par(2))/(par(4)-1)*(0:par(4)-1);
   s1 = -s1;  % fit to TL !
   
   ntrf = size(s1,3); 
   nf   = 1; 
   ifreq= 1;
   nr   = par(7)
   nd   = par(4)
   lss = ones(ntrf,1);
   xval = 45:1:90;
   if ij == 1
      nototal = zeros(nd,nr,length(xval),nfile);
   end
   % the uncertainty map
   for id=1:nd
      for ir=1:nr
         ss = squeeze(s1(id,ir,:));
         [no,xo]= mlpdf(ss,lss,xval);
         nototal(id,ir,:,ij) = no;
      end
   end
   Nototal(:,:,:) = sum(nototal,4);
   
   if iplot & ij == nfile
      figure
      imagesc(ranges,depth,squeeze(s1(:,:,1)),[50 80]),
      colormap(flipud(jet));
      colorbar4
      style1
      
      % test
      pp = nan*ones(nr,length(xval));
      id = 7
      for ir=1:nr
         ss = squeeze(s1(id,ir,:));
         [no,xo] = mlpdf(ss,lss,xval);
         pp(ir,:)= no./(sum(no));
      end

      cbound=[0 0.2];
      figure
      %subplot(2,4,[1:3])
      subplot(3,1,[1 2])
      HH = imagesc(ranges,xval,pp',cbound),
      colormap(clrscl('wk',10))
      axis ij
      %colormap(flipud(gray));
      xlabel('Range (m)'), ylabel('TL (dB)'); 
      axis([min(ranges) max( ranges) 45 85 ])
      titlenew('(a)')
%      style1
      pos_fig = get(gca,'position');
      H = colorbar('horiz','XAxisLocation','top', 'Location', 'NorthOutside');
      pos = get(H,'position');
      set(H,'position',...
            [pos(1)+pos(3)/3*2 pos(2)+pos_fig(4)*0.3 pos(3)/3 pos(4)*0.65])
      text(pos(1)+pos(3)/3*2, pos(2)+pos_fig(4)*0.5,'pdf',...
           'HorizontalAlignment','right','Units','Normalized')
      set(gca,'position',pos_fig)
%      keyboard
 
      ir1 = min(find(ranges>=2100)); 
      ir1 = min(find(ranges>=2750)); 
      hold on; 
      plot(ranges(ir1)*[1 1],[xval(1) xval(end)],'k','linewidth',2)
      plot(ranges(ir1)*[1 1],[xval(1) xval(end)],'w:','linewidth',2)
      ir2 = min(find(ranges>=2300)); 
      ir2 = min(find(ranges>=2850)); 
      hold on; 
      plot(ranges(ir2)*[1 1],[xval(1) xval(end)],'k','linewidth',2)
      plot(ranges(ir2)*[1 1],[xval(1) xval(end)],'w:','linewidth',2)
      
      cbound = [0 0.25];
      subplot(3,2,5)
      sp = spread(pp(ir1,:),xval);
      spreadarr1 = sp(2)-sp(1);    
      med1 = sp(3);lval1 = sp(1);uval1 = sp(2);
      xx1 = sort([xval lval1 uval1]);
      pp1 = interp1(xval,pp(ir1,:),xx1);
      ilmin = find(xx1==lval1);ilmxn = find(xx1==uval1);
      xaxis = xx1([ilmin:ilmxn]);
      fill([lval1 xaxis uval1],[0 pp1([ilmin:ilmxn]) 0],[0.9 0.9 0.9] )
      hold on
      plot(xx1,pp1,'k','linewidth',2)
      plot(med1*[1 1],cbound,'r-')
      ylabel('pdf')
      xlabel('TL (dB)')
      titlenew('(b)')
      %     ylabel('TL (dB)'); 
      axis([45 85 cbound])
      style1
      plot(lval1*[1 1],cbound,'k--')
      plot(uval1*[1 1],cbound,'k--')
      text(lval1*0.98,max(cbound)*0.7,'5-th','HorizontalAlignment',...
               'right');
      text(uval1*1.02,max(cbound)*0.7,'95-th','HorizontalAlignment',...
           'left')

      
      subplot(3,2,6)
      sp = spread(pp(ir2,:),xval);
      spreadarr1 = sp(2)-sp(1);    
      med1 = sp(3);lval1 = sp(1);uval1 = sp(2);
      xx1 = sort([xval lval1 uval1]);
      pp1 = interp1(xval,pp(ir2,:),xx1);
      ilmin = find(xx1==lval1);ilmxn = find(xx1==uval1);
      xaxis = xx1([ilmin:ilmxn]);
      fill([lval1 xaxis uval1],[0 pp1([ilmin:ilmxn]) 0],[0.9 0.9 0.9] )
      hold on
      plot(xx1,pp1,'k','linewidth',2)
      plot(med1*[1 1],cbound,'r-')
      xlabel('TL (dB)')
      titlenew('(c)')
      %     ylabel('TL (dB)'); 
      axis([45 85 0 0.25])
      style1
      plot(lval1*[1 1],cbound,'k--')
      plot(uval1*[1 1],cbound,'k--')
      text(lval1*0.98,max(cbound)*0.7,'5-th','HorizontalAlignment',...
               'right');
      text(uval1*1.02,max(cbound)*0.7,'95-th','HorizontalAlignment',...
           'left')
      print -depsc illustration
   end
end

print -dpng fig5abc

for id=1:nd
   for ir=1:nr
      no = Nototal(id,ir,:);
      sp = spread(no,xo);
      spreadarr(ir,id) = sp(2)-sp(1);    
      med(ir,id) = sp(3);
      lval(ir,id)= sp(1); 
      uval(ir,id)= sp(2);
   end
end
%plot(ranges,uval(:,7))
%%
figure
     subplot(3,1,[1 2])
     idep=7
     fill([ranges fliplr(ranges)],[uval(:,idep)' fliplr(lval(:,idep)')],[0.9 0.9 0.9] )
hold on

plot(ranges,med(:,idep),'r')

plot(ranges,lval(:,idep),'g')
ylabel('TL (dB)'), xlabel('Range (m)')
set(gca,'ydir','rev','xlim',[min(ranges) max(ranges)],'ylim',[40 85])
print -dpng fig5d
