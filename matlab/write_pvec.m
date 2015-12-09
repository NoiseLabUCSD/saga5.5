function status=write_pvec(ppres,freq,depth)
%write_pvec(ppres,freq,depth)
%
% write the pressure vector for  frequency 'freq' to the file 
% pres.in
% pres should 1 column complex vector
%
% Peter Gerstoft
[n,m]= size(ppres);
if (n~=length(depth))
  'Error pressure and depth array do not have consistent dimensions'
  keyboard
end
  
dep=(1:n);
%pg 3 Apr 01: this should not be done here. ppres  =  ppres/(norm(ppres));
fidout =  fopen('pres.in','a');
fprintf(fidout,' pressure vector \n');
fprintf(fidout,' %f\n',freq);
fprintf(fidout,' %d\n',n);
fprintf(fidout,' %f\n',depth);
y= [ dep; real(ppres)';imag(ppres)'];
fprintf(fidout,' %d  (%e,%e) \n',y);
fprintf(fidout,'!  \n');
%status=fclose(fidout)







