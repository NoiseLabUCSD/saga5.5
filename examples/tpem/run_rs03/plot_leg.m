
for iframe=7             %7:20
if (iframe>=7 & iframe<=11)
  gen='./gendata/gendata5';
  load  prof5  % profile 5
elseif   (iframe>=12 & iframe<=13)
  gen='../gendata/gendata6'
  load ../../profiles/prof6  % profile 6
elseif   (iframe>=14 & iframe<=17)
  gen='../gendata/gendata7'
  load  ../../profiles/prof7  % profile 7
elseif   (iframe>=18 & iframe<=20)
  gen='../gendata/gendata8'
  load  ../../profiles/prof8  % profile 8
end
[genrang,genpres,comment] = read_Dformat([ gen '.in']);

if iframe < 10;
   Fname = sprintf('%s%s%s','F040298_0',num2str(iframe),'_2');
else;
   Fname = sprintf('%s%s%s','F040298_' ,num2str(iframe),'_2');
 end;
 inputfil=Fname
  
  plotsaga
%print -djpeg temp

figure(iframe)
clf
subplot(5,1,3) %coverage diagram
[sinv,par]	= readascii([inputfil '.trf']);
jrange =round([1 (par(7)/2) par(7)]);
jreceiv=round([ 1 par(4)*.4 3]);
nplot=3;                           % number of plots for each figure
 height=linspace(par(2),par(3),par(4));
 ranges=par(5)+par(6)*( 0:par(7)-1);
ifr=1;  
[OneFreqinv ] =  coverage2(sinv,par,ifr,[-150 -110]); % range-independent
xlabel('')
%axis([0 100 0 500])
axis('xy'); title(' ')
%text(1,180,['best fit:  ' num2str(bestfit,4)],'fontsize',11)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(5,1,4) %reference coverage diagram
[sgen,par]	= readascii([gen '.trf']);
[OneFreqgen ] =  coverage2(sgen,par,ifr,[-150 -110]); % range-independent
%axis([0 100 0 200])
axis('xy');xlabel('')
title(' ')
%text(1,180,['Reference'],'fontsize',11)
xlabel('Range (km)'); xlabel('')

subplot(5,1,5) %reference coverage diagram
plgen=20*log10(abs(OneFreqgen));
plinv=20*log10(abs(OneFreqinv));
pldiff=abs(plgen-plinv);
medfiltmat=[floor(6/(height(1)- height(2))) floor(10/par(6))]
plmed = medfilt2(pldiff,medfiltmat);
imagesc(  ranges,height,pldiff),h=colorbar;
pos=get(h,'position'); pos(3)=0.025; pos(1)=pos(1)-0.023;
set(h,'position',pos,'xlim',[0.3 1])
imagesc(  ranges,height,plmed,[3 33]),h=colorbar;
pos=get(h,'position'); pos(3)=0.025; pos(1)=pos(1)-0.023;
set(h,'position',pos,'xlim',[0.3 1])
axis('xy');xlabel('')
%text(1,180,['Difference'],'fontsize',11)
xlabel('Range (km)');
ylabel('Height (m)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%subplot(5,1,1)     % inversion match
subplot('Position', [0.13 0.795306 0.65 0.129694])
hold off
plot(a(:,j1+1),a(:,j1+2),'k-',a(:,j1+3),a(:,j1+4),'r','linewidth',1.5)
% for bw plot plot(a(:,j1+1),a(:,j1+2),'b-',a(:,j1+3),a(:,j1+4),'b')
%    axis(plaxis(:,j));
%axis([0 100 -10 100])
ylabel('Clutter (dB)')
%xlabel('range (km)')
hold    on
plot(genrang/1000,genpres+mean(a(:,j1+4))-mean(genpres),'b-','linewidth',1) 
%h=legend('measured','inverted','M-prof  ')
%title(['Spandar F040298, frame' num2str(iframe)])
%set(h,'fontsize',20)
% for bw plot(genrang/1000,genpres+mean(a(:,j1+4))-mean(genpres),'-.') 
%keyboard
%plot(clut(:,1),clut(:,2),'r')
%plot(clutgen(:,1),clutgen(:,2),':g')
%legend('data','best fit','clutter cs', 'clutter cs gen')
axis([ 10 60 0 80])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(5,1,2) % plot profiles 
%eval(['load em_prof' testcase])
nprof=size(ref,1);
for iprof=1:nprof
  ref(iprof,:)=ref(iprof,:)-ref(iprof,2)+330;
end
hold off 
alt=alt*.3048; % now it is in m  
rng=rng*1.852; % now it is in km
pcolor(rng*ones(1,size(alt,2)),alt,ref), caxis([310 350]),shading interp, h=colorbar;
pos=get(h,'position'); pos(3)=0.025; pos(1)=pos(1)-0.023;
set(h,'position',pos,'xlim',[0.3 1])

hold on
    rat = 0.05;	% conversion constant: how many M units per range unit
   for i=1:5:nprof
      plot( rng(i)+rat*(ref(i,:)-(ref(i,1))), alt(i,:),'b-','linewidth',1.5 )
%      plot( rng(i)+rat*(vel(:,2)-(vel(1,2))), vel(:,1),'r','linewidth',1.5 )
    end
ylabel('height (m)')

xx=res(:,1)
thick= xx(8);
mdef = xx(7);
c1   = xx(9);
delta=xx(10); % no evaporation duct
x=-1:0.01:1;
baseh=xx(1);


load markov.in
for ii=2:6
%  baseh=baseh+xx(ii)*x.^(ii-1);
  baseh=baseh+xx(ii)*markov(:,ii)';
  end
baseh=max([baseh; baseh*0]);
xran=(x+1)/2*60;
xran=((markov(:,1)'-1)/99)*60;

plot(xran+rat*c1*baseh,baseh,'r','linewidth',1 )
plot(xran-mdef*rat+rat*c1*baseh,baseh+thick,'r','linewidth',1 )

for ii=1:5
  [minval,I]=min(abs(xran-(-5+12*ii)));
  [mprof,z]=trilin(baseh(I),thick,mdef,c1,delta);
  plot( xran(I)+rat*(mprof-330), z,'r-','linewidth',1.5)
end

axes('position',[0,0,1,1]); axis('off')
bs=0.15; as=.95; da=0.173
text(bs,as,'a)','fontsize',16) 
 text(bs,as-da,'b)','fontsize',16) 
 text(bs,as-2*da,'c)','fontsize',16) 
 text(bs,as-3*da,'d)','fontsize',16) 
 text(bs,as-4*da,'e)','fontsize',16) 


%eval(['print -depsc '  inputfil]) 
%  eval(['print -dpng '  inputfil]) 

end



