% read_mh_bin script
fid = fopen( filename, 'r' );
if fid<0
  display([ 'File ' filename ' does not exist'])
  keyboard
  end
%lrecl = fread( fid, 1, 'long' );

%Try changing this line for the sun, if doesn't work
%lrecl = 4 * lrecl;

%rec = 0;  %Record length is one less than in the krakm.f file
%fseek( fid, rec * lrecl, -1 );
%fseek( fid, 4, 0 );
dummy = fread(fid,1,'int');
nparm = fread(fid,1,'int');
%
dummy = fread(fid,2,'int');
nppd = fread(fid,1,'int');
%rec = 1;
%fseek( fid, rec * lrecl, -1 );
dummy = fread(fid,2,'int');
ndigit= fread(fid,nparm,'int');
%rec = 2;
%fseek( fid, rec * lrecl, -1 );
dummy = fread(fid,2,'float');
f_min = fread(fid,nparm,'float');
%rec = 3;
%fseek( fid, rec * lrecl, -1 );
dummy = fread(fid,2,'float');
f_max = fread(fid,nparm,'float');

%rec = 4;
%fseek( fid, rec * lrecl, -1 );
dummy = fread(fid,2,'float');
df = fread(fid,nparm,'float');
%
dummy = fread(fid,2,'float');
par2phy = fread(fid,nparm,'int');
%
dummy = fread(fid,2,'float');
iopt = fread(fid,40,'int');

nobs=prod(ndigit);   
 nparmp1= nparm+1;
%%%%% xtt = zeros( nobs, nparmp1 );   %number of modes
%keyboard
% for ii = 1: nobs
%%    rec = 4 + ii;
%%    fseek( fid, rec * lrecl, -1 );
%    dummy = fread(fid,2,'float');
%    xdum = fread( fid, nparmp1 , 'float' )'; %Data is read columwise first
%    xtt(ii,:)=xdum;
%  end
xtt = fread( fid, [nparmp1+2 inf] , 'float' ); % The first two columns


fval=xtt(end,1:end-1)'; %;	fval=fval/sum(fval);

xtt=xtt(3:(end-1),1 :end-1)';
fclose( fid );
