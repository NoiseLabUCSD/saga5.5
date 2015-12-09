function status=write_covmat(freqvec,Cov,depth)
%write_covmat(freqvec,Cov,depth)
%
% Write Covariance Matrices to file cov_dpss.in
%
% freqvec = vector of frequencies in Hertz
% Cov     = matrix of covariance matrices
%           this must be of size [n,n, length(freqvec)]
% depth   = depth of sensors, length(depth)=n
[n,m,F]= size(Cov);
if(m ~= n)
   fprintf(1,'size(Cov) = [%d %d], length(freqvec)=%d\n', ...
           size(Cov),length(freqvec));
   error('sizes of Cov and freqvec are not compatible');
end
if(F ~=length(freqvec))
   fprintf(1,'size(Cov) = [%d %d], length(freqvec)=%d\n', ...
           size(Cov),length(freqvec));
   error('sizes of 3Cov and freqvec are not compatible');
end

if(n ~= length(depth))
   fprintf(1,'size(Cov) = [%d %d], length(freqvec)=%d\n', ...
           size(Cov),length(freqvec));
   error('sizes of Cov and depth are not compatible');
end
%depth=linspace(18.72,112.72,n);

fidout =  fopen('cov.in','a');

for ifreq=1:F
   fprintf(fidout,' estimated covariance matrices using dpss\n');
   fprintf(fidout,' %f     0.000 dB\n',freqvec(ifreq));
   fprintf(fidout,' %d\n',n);
   fprintf(fidout,' %f\n',depth);
%   cols = [ (ifreq-1)*n+1:ifreq*n ];
   Cx = Cov(:,:,ifreq); % cut one cov-matrix out of the bunch
   for row=1:n
     fprintf(2,'%3d\b\b\b',row);
      for col=1:n
         fprintf(fidout,'%10d%10d (%E,%E) \n',row,col, ...
                 real(Cx(row,col)),imag(Cx(row,col)));
      end
   end
end
fprintf(fidout,'!  \n');
status=fclose(fidout);
unix('ls -l cov.in');



