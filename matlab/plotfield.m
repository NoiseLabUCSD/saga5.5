function [st,pulse,t,depth]=plotfield(filename,strng,frq,dpth_rng,suboption,col);
% [st,pulse,t,depth]=plotfield(filename,strng,frq,dpth_rng,suboption,col);
%     strng = tlran, tldep, timed, timer
%     suboption env
%----------------------
% Plot output generated
% by TPEMBB
% Peter Gerstoft Oct 96/ june 00
%----------------------
% read the output from tpem/saga/oases/prosim, written to file "filename"
% s contains the response for each range depth and frequency. They are
% all packed into the first index of s.
% par contains the parameters used in the processing
% par(1) = 1 
% par(2) = first receiver height
% par(3) = last receiver height
% par(4) = number of receivers in height
% par(5) = first range
% par(6) = delta range
% par(7) = number of ranges
% par(8) = nx, number of time samples, T=nx*(delta t)
% par(9) = lx, first sample computed by tpem
% par(10) = mx, last sample computed by tpem
% par(11) = delta t
% par(12) = 1 
%
clear readascii;
st=[]; pulse=[]; t=[]; depth=[];
[s,par]	= readascii(filename);

%s=exp(-1i*3/4)*s;
                                  % based on the tpem run here is
fs=1/par(11);			  % sampling frequency=1/delta t
df=1/(par(11)*par(8));            % delta freq=1/tmax
%par(11)
%par(8)
fmin=df*(par(9)-1);               % the minimum used frequency
fmax=df*(par(10)-1);              % the maximum used frequency
nrec=par(4);
nrng=par(7);
%                                   %height range
range=par(5)+par(6)*(0:nrng-1);
del_range=par(6);
if(nrec>1)
  del_depth=(par(3)-par(2))/(nrec-1);
  depth=(par(3)-par(2))/(nrec-1)*(0:nrec-1)+par(2);
else
  del_depth=max([10 par(3)-par(2) par(2)]);
  depth=par(3);
end;
%
t=0:par(11):(par(8)-1)*par(11);
freq=0:df:df*(par(8)/2);
nf=par(10)-par(9)+1;              % number of used fequencies
ifrq=round(frq-fmin)/df+1;

tp=.100;

dpth_rng
han=myhanning(fmin,fmax,freq,df);
if(strng == 'tlran')
   depth
   idpth=find(depth == dpth_rng);
  idpt=find(depth < dpth_rng+del_depth/2 & depth > dpth_rng-del_depth/2)
   for ii=1:length(ifrq)
     plot(range,20*log10(abs(s(ifrq(ii),idpth:nrec:nrng*nrec))),col)
     hold on;
     plot(range,20*log10(abs(s(ifrq(ii),idpth:nrec:nrng*nrec))),'r*')
   end
   hold off
elseif(strng == 'tldep')
   irng=find(range == dpth_rng)
  if (size(irng,1)==0) 
    disp ' no matching range'
    disp 'available ranges are' 
	  range
    disp 'requested interval'
	  [dpth_rng-del_range/2  dpth_rng+del_range/2]
	  return
  end
   for ii=1:length(ifrq)
     plot(depth,20*log10(abs(s(ifrq(ii),[(irng-1)*nrec+1:irng*nrec]))),col)
     hold on;
     plot(depth,20*log10(abs(s(ifrq(ii),[(irng-1)*nrec+1:irng*nrec]))),'r*')
   end
   hold off
elseif(strng=='timer')
  idpt=find(depth < dpth_rng+del_depth/2 & depth > dpth_rng-del_depth/2)
  if (size(idpt,1)==0) 
    disp ' >>>>> no matching depth. Available depths are' 
	  depth
    disp 'Based on input depth must be in the interval' 
	  [dpth_rng-del_depth/2 dpth_rng+del_depth/2]
	  return
  end
%keyboard
  eval(['ss([1:nrng],1:par(8)/2+1)=' ...
        '[zeros(nrng,par(9)-1) s(:,(idpt-1)*nrng+[1:nrng])''' ...
        'zeros(nrng,par(8)/2-par(10)+1)];']);
  wtap=costap(par(8)/2,par(9),par(10),2);
%  rwx=rxx.*wtap;
%  rwx=[wtap 0];
  rwx=han;
  rxx=rwx;
%  rxx=[rwx 0. conj(rwx(length(rwx):-1:2))];
%  rxx=[rwx 0.];
  sourxx=real((ifft([rxx conj(rxx(length(rxx)-1:-1:2))])));
  ss1=ss.*(ones(nrng,1)*rxx(:).');
  ss2=[ss1 conj(ss1(:,size(ss,2)-1:-1:2))];
%  sstmp1=ss(3,:).*rxx;
%  sstmp2=[sstmp1 conj(sstmp1(length(sstmp1)-1:-1:2))];
  for ii=1:nrng
    st(ii,:)=real(ifft(ss2(ii,:)));
    if(suboption=='env')
      st(ii,:)=abs(hilbert(st(ii,:)));
    end
  end;
  tim=[0:par(11):(par(8)-1)*par(11)];
  t=tim;
  depth=range;
  pulse=[sourxx];
%keyboard
%  hwigglexy3(st,[1/32768 0.02],[0 .02],15,0,tim,range,'r-');
%  hwiggle(st,[1/32768 20],[0 20],2,0);
elseif(strng=='timed')
  irng=find(min(range-dpth_rng))
  irng=find(range < dpth_rng+del_range/2 & range > dpth_rng-del_range/2)
  if (size(irng,1)==0) 
     disp ' no matching range'
    disp 'available ranges are' 
	  range
    disp 'range must be in the interval' 
	  [dpth_rng+del_range/2  dpth_rng-del_range/2]
	  return
  end
%  for ii=1:nrec
    eval(['ss([1:nrec],1:par(8)/2+1)=' ...
        '[zeros(nrec,par(9)-1) s(:,(irng-1)*nrec+[1:nrec])''' ...
        'zeros(nrec,par(8)/2-par(10)+1)];']);
%  end;
%  [chirp_time,chirp_dum,chirp_spec]=chirp(par(8),fs,tp,fmin,fmax);
%  [autospec]=conj(chirp_dum).*(chirp_dum);
%  rxx=autospec(1:length(autospec)/2+1).*han;
%  [rxx,txx]=suscharge(0.82,1,90.0,fs,par(8),4);
%  enrgyt=sum(txx.^2)/fs
%  enrgyf=sum(2*(abs(rxx)/fs).^2)*fs/par(8)
%  wtap=costap(par(8)/2,par(9),par(10),2);
%  rwx=rxx.*wtap;
%  rwx=wtap;
%  rxx=[rwx 0. conj(rwx(length(rwx):-1:2))];
%  rxx=[rwx 0.];
  rxx=han;
  sourxx=real((ifft([rxx rxx(length(rxx)-1:-1:2)])));
  pulse=[sourxx];
  ss1=ss.*(ones(nrec,1)*rxx(:)');
  ss2=[ss1 conj(ss1(:,size(ss,2)-1:-1:2))];
%  sstmp1=ss(3,:).*rxx;
%  sstmp2=[sstmp1 conj(sstmp1(length(sstmp1)-1:-1:2))];
  for ii=1:nrec
    st(ii,:)=real(ifft(ss2(ii,:)));
    if(suboption=='env')
      st(ii,:)=abs(hilbert(st(ii,:)));
    end
  end;
  t=[0:par(11):(par(8)-1)*par(11)];
%  hwigglexy3(-hst,[1/fs 1],[0 0],1.2,0,t*1000,depth,'r-')
%  set(gca,'ydir','reverse')
end;






