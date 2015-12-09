%SAGA_IN_GEN
%       
%  File to generate input file for saga. This produces a file called test.in
% which conforms to the required input structure to run saga programs.
%
% USAGE
%
%  >> saga_in_gen(FREQ,depth,prr);
%
%   Where the following variables are:
% 
%   FREQ       ==> Vector containing the frequecies
%   depth      ==> Vector containing hydrophone depths
%   prr        ==> Pressure field matrix (rows are depth and columns are 
%                  frequencies
%

function saga_in_gen(FREQ,depth,prr,opt)

fid = fopen('test.in','w');
fprintf(fid,'%c', '! File format: Hydrophone vectors');
fprintf(fid,'\n');
fprintf(fid,'%c', 'Vertical array data');
fprintf(fid,'\n');
for mm = 1:length(FREQ)
 fprintf(fid,'%4.3f %i %i\n', [FREQ(mm) length(FREQ) 1]);
 fprintf(fid,'%i\n', length(depth));
 for nn = 1:length(depth)
  fprintf(fid,'%4.2f\n', depth(nn));
 end;
 for nn = 1:length(depth)
  if (opt == 'magonly')
  strng = [num2str(nn),' (',num2str(abs(prr(nn,mm)),15),', ',num2str(0),'), ',num2str(abs(prr(nn,mm)),15)];
  else
  strng = [num2str(nn),' (',num2str(real(prr(nn,mm)),15),', ',num2str(imag(prr(nn,mm)),15),'), ',num2str(abs(prr(nn,mm)),15)];
  end;
  fprintf(fid,'%c', strng);
  fprintf(fid,'\n');
 end;
  fprintf(fid,'%c', ['Next set of data']);
  fprintf(fid,'\n');  
end;
 








