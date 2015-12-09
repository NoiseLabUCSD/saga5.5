function  nfile = writefort81(filenm,mm)
% generate fort.81 from MCMC samples
%filenm = 'asiaexcsdm135';
%   mm = 10000;

filename=[filenm '.mh1'];
read_mh_bin;
fitval = fval;
names = dir(filename);
if names.bytes < 90E6;
   filename = [filenm '.mh2'];
   read_mh_bin;
   fitval = [fitval; fval];
 else
   display('Warning *mh2 file not read!')
   keyboard
end

[fitval ind] = sort(fval);
mparm = size(xtt,2);       % number of the parameters
nsample = length(fitval) ; % number of m's
nfile = floor(nsample/mm); % number of fiels
nfile = 1;
jj = rand(nsample,1);
jj = floor(jj./max(jj)*nsample);
if find(jj == 0) 
   jj(find(jj == 0)) = 1;
end
samp = xtt(ind(jj),:); xtt = [];
fitval = fval(ind(jj));fval = [];

eval([filenm '_mh']);
iforward = iopt(30);
% move geometrical parameters to the end: source range, array position
orderparm = [1: mparm];
%keyboard
if ~isempty(find(par2phy == 9 | par2phy== 8 | par2phy == 21) )
   endind = [];   
   if ~isempty(find(par2phy == 9))
      endind = [endind find(par2phy == 9).'];   
   end
   if ~isempty(find(par2phy == 8))
      endind = [endind find(par2phy == 8).'];   
   end
   if ~isempty(find(par2phy == 21))
      endind = [endind find(par2phy == 21).'];   
   end
   for ii = 1:length(endind)
      orgind = orderparm(find(orderparm ~= endind(ii))); 
      orderparm = [orgind endind(ii)];      
   end
end

for ii = 1:nfile
   fidout = fopen(['fort.81_' num2str(ii)],'w');
   for jj = (ii-1)*mm+[1:mm];
      fprintf(fidout,' %f',[samp(jj,orderparm) fitval(jj)] );
      fprintf(fidout,' \n');
   end 
   status = fclose(fidout);
end

