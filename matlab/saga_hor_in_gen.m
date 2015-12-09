%SAGA_HOR_IN_GEN
%       
%  File to generate input file for saga. This produces a file called test.in
% which conforms to the required input structure to run saga programs.
%
% USAGE
%
%  >> saga_hor_in_gen(FREQ,depth,range,prr,opt);
%
%   Where the following variables are:
% 
%   FREQ       ==> Vector containing the frequecies
%   range      ==> Vector containing hydrophone ranges
%   prr        ==> Pressure field matrix (rows are  and columns are 
%                  frequencies
%

function saga_in_gen(FREQ,depth,range,prr,opt)
%keyboard
fid = fopen('test.in','w');
fprintf(fid,'%c', '! File format: Hydrophone vectors');
fprintf(fid,'\n');
fprintf(fid,'%c', 'Horizontal array data');
fprintf(fid,'\n');
for nn = 1:length(range)
 for mm = 1:length(FREQ)
 fprintf(fid,'%4.3f %i %i\n', [FREQ(mm) length(FREQ) nn]);
 fprintf(fid,'%i\n', length(depth));
 for oo = 1:length(depth)
  fprintf(fid,'%4.2f\n', depth(oo));
 end;
  if (opt == 'magonly')
  strng = [num2str(1),' (',num2str(abs(prr(nn,mm)),15),', ',num2str(0),'), ',num2str(abs(prr(nn,mm)),15)];
  else
  strng = [num2str(1),' (',num2str(real(prr(nn,mm)),15),', ',num2str(imag(prr(nn,mm)),15),'), ',num2str(abs(prr(nn,mm)),15)];
  end;
  fprintf(fid,'%c', strng);
  fprintf(fid,'\n');
  fprintf(fid,'%c', [' ! Comments']);
  fprintf(fid,'\n');  
 end;
%  fprintf(fid,'%c', ['Next set of data']);
%  fprintf(fid,'\n');  
end;
 








