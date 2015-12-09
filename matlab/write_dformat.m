function status=write_dformat(pres,comment)
%write_Dformat(pres,comment)
%%%%%
% pres(frequecies, range,depth)
% write the pressure vector for  ranges 'range' to the file 
% pres.in
%
% Peter Gerstoft
[nfr,nr,nd]= size(pres);

fidout =  fopen('pres.in','w');
fprintf(fidout,'! %s\n',comment);
fprintf(fidout,'! frequencies range depth pressure \n','');
for ifr=1:nfr
  for ir=1:nr
    for id=1:nd
      
      fprintf(fidout,' %d %d %d   (%e, %e) \n',...
      [ifr id ir real(pres(ifr,ir,id)) imag(pres(ifr,ir,id))]);
    end
  end 
end
%keyboard
fprintf(fidout,'!  \n');
%status=fclose(fidout)







