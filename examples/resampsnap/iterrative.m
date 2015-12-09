filenm   ='asiaexcsdm135.trf'

[s1,par] = rereadtrf(filenm);
ranges   = par(5)+par(6)*( 0:par(7)-1);
depth    = par(2)+ (par(3)-par(2))/(par(4)-1)*(0:par(4)-1);
% don't run the following line. Memory problem
%s_tl     = -20*log10(abs(s1));  % fit to TL !

ntrf = size(s1,3);
nf   = 1; ifreq=1;
nd   = par(4)
nr   = par(7) 


if (1==-1)
figure
subplot(2,1,1)
imagesc( ranges,depth,squeeze(-s_tl(:,:,1))),colorbar;
subplot(2,1,2)
imagesc( ranges,depth,squeeze(s_tl(:,:,ntrf))),colorbar;
colormap(fliplr(jet))
data = squeeze(s_tl(:,:,ntrf));
figure;imagesc(data);
end

figure
[xval,marg,fitval] = plotenum('asiaexcsdm135',3);
xlabel('Sediment sound speed (m/s)')
expfit = exp(-fitval);
expfit = expfit(1:ntrf);
lss    = expfit/sum(expfit);
%print -depsc marg2d

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% here we check data, 
% problem 1: SD is not correct (20m?)
%%%%%%       That value is for the tl uncertainty
%
% problem 2: RD is not right (99.5 24.5 15), It's not evenly distributed.
%            Ele#4 is dead. Further, the array tilt 
% problem 3: I intented to regenerate the replica using peter's resamp.f
%            but the tl is not right ???
% for dead ele (#4), maybe I should fill in zero at the element in d. 
%            done. The data are in the following directoty. infor and
%            csdm*** 
%            in dat file, rd: 99.5 24.5 16 tilt 
datapath = '~chenfen/ASIAEX/matlab.dir/fftfx/csdm_mod_16.dir/';
eval(['load ',datapath,'infor']);
fff      = [95 195 295 395 805 850 905];       % signal freguencies 
ifre     = 3;  % select 295 Hz;
nsnap    = 4;
deptrial = linspace(45,55,20);
nenv     = 40*40;
Neff     = 5;
   %--- remind myself replica structure (s1)
   [mele, mrng, menv] = size(s1);
   BB = zeros(mrng,menv);
%
 for irng = 1:mrng
      for ienv = 1:menv
         repl = s1(:,irng,ienv);
        s1(:,irng,ienv)= s1(:,irng,ienv)/norm(repl);
      end
   end

%--- load csdm
datablocks = 95:150  
%datablocks = [130 135 140 145];  
ndata=length(datablocks)
irange=round((sr_save(datablocks(1))-ranges(1))/par(6))
[xx idepth]=min(abs(sd_save(datablocks(1))-deptrial))
% main processing loop
comploop

figure
for idata = 1:ndata
%  subplot(4,4,idata)
  subplot(7,8,idata)
  imagesc(xval(3,:),xval(2,:),reshape(obj(idata,:),size(xval,2),size(xval,2) ));
  if (rem(idata,8)~=1),   set(gca,'yticklabel',[]);  end
  set(gca,'xtick',[min(xval(3,:)) max(xval(3,:))])
  if (idata<8*6+1),    set(gca,'xticklabel',[]);  end
end
subplot(7,8,1)
title('Single dataframes')
% print -depsc single

figure
obj1=ones(size(obj(idata,:)));
for idata = 1:ndata
  obj1=obj1.*obj(idata,:);
  subplot(7,8,idata)
  imagesc(xval(3,:),xval(2,:),reshape(obj1,size(xval,2),size(xval,2) ));
  if (rem(idata,8)~=1),   set(gca,'yticklabel',[]);  end
  set(gca,'xtick',[min(xval(3,:)) max(xval(3,:))])
  if (idata<8*6+1),    set(gca,'xticklabel',[]);  end
end
subplot(7,8,1)
title('Forward likelihood')
% print -depsc forward

 
%figure
%obj1=ones(size(obj(idata,:)));
%for idata = 1:ndata
%  obj1=obj1.*obj(idata,:);
%  if (idata>5),  obj1=obj1./obj(idata-5,:); end
%  subplot(7,8,idata)
%  imagesc(xval(3,:),xval(2,:),reshape(obj1,size(xval,2),size(xval,2) ));
%end
%subplot(7,8,1)
%title('Forward likelihood, max 5 data blocks')
% print -depsc forward5

 

%%%% more testing

figure(2); clg;
figure(3); clg;
%%%%%% new likelihood
phi1=zeros(size(phi(1,:)));
for idata = 1:ndata
  phi1=phi1+ phi(idata,:);
  nu=min(min(phi1))/Neff
  obj1=exp(-phi1./nu);
  obj1=obj1.*obj1;  % neff
  obj1=obj1./sum(sum(obj1));
 figure(2); subplot(7,8,idata)
  imagesc(xval(3,:),xval(2,:),reshape(obj1,size(xval,2),size(xval,2) ))%,colorbar;
  if (rem(idata,8)~=1),   set(gca,'yticklabel',[]);  end
  set(gca,'xtick',[min(xval(3,:)) max(xval(3,:))])
  if (idata<8*6+1),    set(gca,'xticklabel',[]);  end

  lh=reshape(obj1,size(xval,2),size(xval,2) );
  [l1 i1]=max(lh);[l2 i2]=max(l1);
  mlval(idata,1)=xval(3,i2); mlval(idata,2)=xval(2,i1(i2));
  likeh1(idata,:,:)= lh;
  sp1=spread( lh(i1(i2),:),xval(3,:)); sp2=spread( lh(:,i2),xval(2,:));
  lb(idata,:)=[sp1(1) sp2(1)];       ub(idata,:)=[sp1(2) sp2(2)];
figure(3)
subplot(1,2,1)
hold on
maxlike(idata,:)=lh(i1(i2),:); plot(xval(3,:),lh(i1(i2),:))
xlabel('Sound speed (m/s)')
subplot(1,2,2)
hold on
maxlike2(idata,:)=lh(:,i2)'; plot(xval(2,:),lh(:,i2))
xlabel('Bathymetry (m)')
  end
% print -depsc iterative

figure(3), style1
print -deps maxline
cval=[0 5e-3];
figure; subplot(1,2,1) 
imagesc(xval(3,:),1:56,maxlike,cval); colorbar
ylabel('data block #'); xlabel('Sound speed (m/s)')
style1
subplot(1,2,2) 
imagesc(xval(2,:),1:56,maxlike2,cval); colorbar
ylabel('data block #'); xlabel('Bathymetry (m)')
style1
print -depsc maxcont

% source and range
figure
subplot(2,1,1)
plot(Time(datablocks)/60,sd,'b--',Time(datablocks)/60,sd_save(datablocks),'b-')
ylabel('Source depth (m)');
xlabel('Time (s)')
subplot(2,1,2)
plot(Time(datablocks)/60,sr,'b--',Time(datablocks)/60,sr_save(datablocks),'b-')
ylabel('Source range (m)');
xlabel('Time (s)')
% print -deps location



%%% parameter tracking

for idata = 1:ndata
  lh=reshape(obj(idata,:),size(xval,2),size(xval,2) );
  [l1 i1]=max(lh);[l2 i2]=max(l1);
  mlval(idata,1)=xval(3,i2); mlval(idata,2)=xval(2,i1(i2));
  likeh1(idata,:,:)= lh;
  sp1=spread( lh(i1(i2),:),xval(3,:)); sp2=spread( lh(:,i2),xval(2,:));
  lb(idata,:)=[sp1(1) sp2(1)];       ub(idata,:)=[sp1(2) sp2(2)];
%  imagesc(xval(3,:),xval(2,:),lh)
%  title(num2str(tid(idata)))
%  pause
end


tid=Time(datablocks)/60;


tid=Time(datablocks)/60;

figure
subplot(2,1,1)
fill([tid fliplr(tid)],[ub(:,1)' fliplr(lb(:,1)')],[0.9 0.9 0.9])
hold on;
plot(tid,mlval(:,1),'linewidth',2)
ylabel(' Sediment sound speed (m/s)','fontsize',12)
set(gca,'ylim',[1550 1750])

subplot(2,1,2)
fill([tid fliplr(tid)],[ub(:,2)' fliplr(lb(:,2)')],[0.9 0.9 0.9])
hold on
plot(tid,mlval(:,2),'linewidth',2)
ylabel(' Bathymetry (m)','fontsize',12)
xlabel('Time (s)','fontsize',12)
set(gca,'ylim',[100 120])
% print -deps like5val
% print -deps likeiter
% print -deps likesingle