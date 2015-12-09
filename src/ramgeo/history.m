casedir='/heart5/gerstoft/itworkshop/tc0/'
[pr par]  =read_it([casedir 'ivwkt0_v_0025.cpr']);

[prh parh]=read_it([casedir 'ivwkt0_h_0025.cpr']);

range=parh(2)+parh(3)*[0:(parh(4)-1)];
depth=par(5)+par(6)*[0:(par(7)-1)];
%%%%
figure  % depth 500 m range
subplot(1,2,1)
plot(20*log10(abs(pr(:,1))),depth)
hold on
plot(20*log10(abs(z(20:80,100))),depth,'r')
xlabel('TL (dB)'); ylabel('Depth (m)')
axis([-57 -48 20 80])
title('Vertcal array 25Hz, range 500 m')
subplot(1,2,2)
ff=exp(i*25*2*pi/1488.6*range(100));
plot(unwrap(angle(pr(:,1))),depth)
hold on
plot(unwrap(angle(z(20:80,100).'*ff)),depth,'r')
xlabel('angle (rad)'); ylabel('Depth (m)')
axis([-5 -1 20 80])
norm=(z(20:80,100)*ff)'*(z(20:80,100)*ff)*(pr(:,1)'*pr(:,1))
bartlett=abs((z(20:80,100)*ff)'*pr(:,1))^2/norm
figure 
subplot(2,1,1)
plot(range,20*log10(abs(prh(1,:))))
hold on
plot(range(1:600),20*log10(abs(z(25,1:600))),'r')
ylabel('TL (dB)'); xlabel('Range (m)')
axis([0 3000 -80 -20 ])
subplot(2,1,2)
title('Horzontal array 25 Hz, 25 m depth')
ff=exp(i*25*2*pi/1488.6*range(1:600));
plot(range(1:600),angle(prh(1,1:600)./ff))
hold on
%plot(range(1:100),angle(pram(25,1:100).*ff),'r')
plot(range(1:600),angle(z(25,1:600)),'r')
axis([0 3000 -3.2 3.2])
ylabel('Angle - k_0r (Rad)'); xlabel('Range (m)')

%%%%%


%%%%%%% ramgeo
fid=fopen('tlgeo.line','r');
%lx=fread(fid,1,'int')
tllinegeo2=fscanf(fid,'%f %f',[2 100]);
size(tllinegeo2)
fclose(fid);
plot(tllinegeo2(1,:),-tllinegeo2(2,:),'g')

%%%% reading fort99
fid1=fopen('fort.99','r')
%temp=fscanf(fid1,'%d %f ( %f, %f )',[4 199*600]);
temp=fscanf(fid1,'%d %f  %f %f',[4 199*600]);
pram=temp(3,:)+i*temp(4,:);
pram=reshape(pram,[199 600]);   %pram(depth,range)
tlram=20*log10(abs(pram));
plot(range(1:600),tlram(25,:),'m')
fclose(fid1)
%%%%
disp('reading the ambiguity file');
fid=fopen('../ram_hc/tl.grid','r');
x=fread(fid,[2*m n],'float'); %size(x)
fclose(fid); 
x([1 2*m],:)=[]; x=reshape(x,[2*(m-1) n]);  % delete first and last
                                            % number read by fread
 zhc=x(1:2:(2*m-3),:) + i*x(2:2:(2*m-2),:); % complex number

pramcopy=pram;
tlramcopy=tlram;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prh(1,100)  - pr(6,1) % 500 m range
prh(1,200)  - pr(6,2) % 1000 m range
prh(1,300)  - pr(6,3) % 1500 m range
prh(1,400)  - pr(6,4) % 2000 m range
prh(1,500)  - pr(6,5) % 2500 m range
prh(1,600)  - pr(6,6) % 3000 m range

basename='ivwkt0_h_0'
freq=25:199
pfreq=[]
df=1
for f=freq
  if (f<100)
    filename=[basename '0' num2str(f) '.cpr'];
  else
    filename=[basename num2str(f) '.cpr'];
  end  
  [prh parh]=read_it(filename);
  pfreq=[pfreq; prh(1,:)];
end
dsize=512*3
 tmp = zeros(dsize/2,600);
 tmp(freq*3,:) = pfreq;
 tmp = ([tmp; zeros(1,600); conj(tmp(dsize/2:-1:2,:))]);

 time=linspace(0,3,512*3)
 dt=1/512;
ptime=real(fft( tmp));
range=parh(2)+parh(3)*[0:(parh(4)-1)];

hwiggle(ptime(:,1:10:600)',[dt 50 ],[0 5], 100 , 0,'b')

tred=0.8 +1/1500*range
[ampmax,tmax,pulse,t]=maxamp(ptime(:,1:100),time,tred,[-.5 .5]);
hwiggle(pulse(:,10:5:100)',[dt/10 50 ],[0 5], 100 , 0,'b')

[ampmax,tmax,pulse,t]=maxamp(ptime(:,1),time,1,[-.5 .5]);
hwiggle(pulse(:,100:10:600)',[dt/10 50 ],[0 5], 100 , 0,'b')
axis([0 .1 0 3000])
