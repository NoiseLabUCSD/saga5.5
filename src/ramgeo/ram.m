
Freq=25; 		% source frequency (Hz)

% geometry parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

width=3000; deltar=5;  	        % range in ambiguity space (m)
d1=20; nsd=length(d1);	% VRA element depth (m)
d=1:(200-1);  nd=length(d);  % PE RAM ourput depths
vra_d=114; dep=114;			% water depth at SRA (m) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sed=vra_d+2.5;				% depth at sediment (m)
att1=0.13; att2=0.15;			% attenuation in sediment layer (dB/l
ranges=[0]	% nominal source range
 
rmax=ranges+width;	% maximum range
r=deltar:deltar:rmax; ir=(ranges-width)/deltar; nrr=length(r);
dr=deltar/20; nr=nrr-ir+1;

z=zeros(nd,nr);		% initialize the 3-D matrix

is=1		% VRA depth

m=length(d)+1; n=nrr; 	% size of data

% read the ambiguity file

disp('reading the ambiguity file');
fid=fopen('tl.grid','r');
x=fread(fid,[2*m n],'float'); %size(x)
%x=fread(fid,[2*m n],'double'); %size(x)
fclose(fid); 
x([1 2*m],:)=[]; x=reshape(x,[2*(m-1) n]);  % delete first and last
                                            % number read by fread
 z=x(1:2:(2*m-3),:) + i*x(2:2:(2*m-2),:); % complex number

figure
imagesc(20*log10(abs(z)))

fid=fopen('tl.line','r');
%lx=fread(fid,1,'int')
tlline=fscanf(fid,'%f %f',[2 n]);
size(tlline)
fclose(fid); %x=reshape(x,[2*m n]);
figure
plot(tlline(1,:),-tlline(2,:),'g')


