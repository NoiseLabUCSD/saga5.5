function covdiag2(freqbin,frequency,Cov_Matrix,z)
%diagplot(freqbin,frequency,Cov_Matrix,depth)
%  
%  pretty-print plot of characteristics of a covariance matrix
%  used by plotcov for debugging the signal processing

% swap the variables 'freqbin' and 'frequency' if obviously 
% entered in wrong order:
clf
n=length(Cov_Matrix);
if(freqbin ~= round(freqbin))
  aux = freqbin; freqbin = frequency; frequency = freqbin;
end

displayinfo=[' Frequency' num2str(frequency,6),' Hz, Bin ',num2str(freqbin,6),''];
%uicontrol('Style','text','Position',[490 145 150 20],'String',displayinfo,'Bac%kgroundcolor',[0 1 1]);

subplot(2,3,3)
%title(displayinfo)

     [U,SVD,V]=svd(Cov_Matrix);SVD=diag(SVD);
     ndep=size( SVD,1);
     noise=(sum(SVD(2:ndep)))/(ndep-1)
     SNR=10*log10((SVD(1)-noise)/noise)

els=2; c=1500; omega=2*pi*frequency;
theta = (-90:.25:90);
proj = els * sin(theta .* pi/180.0) .* omega/c;
r = exp(i*(1:n)'*proj)/sqrt(n);
r=r.*(kaiser(n,1.5*pi)*ones(1,size(r,2)));
for ib=1:length(theta)
  beam(ib)=real(r(:,ib)'*Cov_Matrix*r(:,ib));
end
beam=beam/max(beam); beam2=real(r'*U(:,1)/max(r'*U(:,1))).^2;
plot(theta,10*log10(abs(beam)),theta,10*log10(abs(beam2)))
legend('Cov mat','1st EV')
%keyboard
set(gca,'xlim',[-90 90]);
set(gca,'ylim',[-60 0]);
xlabel('Angle (deg)')
%text(-90,1.1,['Freq ' num2str(frequency,4),''])
%text(-900,1.0,['     Bin ',num2str(freqbin,2),''])
%r

title(['Freq ' num2str(frequency),' '])


subplot(2,3,4);
%     SVD_mat=10*log10(svd(Cov_Matrix))+13.98+138+20;     
     SVD_mat=10*log10(svd(Cov_Matrix));
     plot(SVD_mat,'ob');
     eig=SVD_mat(1,1);
%    eig=SVD_mat(1,1)+13.98+138-20*log10(7);
     s = sprintf(' Estimated SNR %6.2f ',SNR);
     title(s,'fontsize',11);
     ylabel('Power (dB)','fontsize',10);
     xlabel('Order of eigenvalues','fontsize',10);
     set(gca,'fontsize',8);
subplot(2,3,1);
     [nrow ncol]=size(Cov_Matrix);
     data=abs(Cov_Matrix);
     data=[data data(:,ncol)];
     data=[data;data(nrow,:)];
     pcolor(data); title('Abs(Cov)','fontsize',12);
     xlabel('hydrophone #','fontsize',10);
     ylabel('hydrophone #','fontsize',10);	
     shading flat
     axis('ij');
     set(gca,'fontsize',8);     
%     colortes
subplot(2,3,2); 
     data=angle(Cov_Matrix);
     data=[data data(:,ncol)];
     data=[data;data(nrow,:)];
     pcolor(data*180/pi); title('Phase(Cov)','fontsize',12);
     xlabel('hydrophone #','fontsize',10);
     ylabel('hydrophone #','fontsize',10);	
     shading flat
     axis('ij');
%     colortes
     set(gca,'fontsize',8);    
subplot(2,3,5);
     [U,S,V] = svd(Cov_Matrix);

     [nsize ndum]=size(Cov_Matrix);
     zvec = linspace(1,nsize,nsize)';
%	z=[9.15 12.15 15.15 18.15 21.15 24.15 27.4 42.4 45.4]; 	
%     z=zvec
     pvec = V(:,1);
     plot(abs(pvec),z,'b*',abs(pvec),z,'b');
%     axis([0 1 0 50]);
     set(gca,'ydir','reverse');
     title('e) Abs EV','fontsize',12);
     ylabel('depth (m)','fontsize',10);
     xlabel('nor. amplitude','fontsize',10);
     set(gca,'fontsize',8);    
%     grid on;
subplot(2,3,6);
     theta=unwrap(angle(pvec))*180/pi;
     plot(theta,z,'m*',theta,z,'m');
%     set(gca,'ylim',[0 50]);
     set(gca,'ydir','reverse');
     title('f) Phase EV','fontsize',12);
     xlabel('angle (deg)','fontsize',10);
     ylabel('depth (m)','fontsize',10);
     set(gca,'fontsize',8);    
%     grid on;
return;






