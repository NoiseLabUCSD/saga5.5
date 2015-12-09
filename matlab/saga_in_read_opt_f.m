function [saga_depth,saga_pr] =  saga_in_read(fname)

if(fname == [])
  fname = input('Input filename ? (.obs file) ','s');
end;
x = findstr(fname,'.obs');
if (x == [])
 fname = [fname,'.obs'];
end;
fid = fopen(fname,'r');
if (fid == -1)
 disp('************No file found***********')
 break
end;
line = fgetl(fid);
line = fgetl(fid);
strng = setstr(fgetl(fid));
eval(['x = [',strng,'];']);
FREQ = [ x(1)];
nfreq = x(2);
nrange = x(3);
strng = setstr(fgetl(fid));
eval(['x = [',strng,'];']);
nphone = x;
for mm = 1:nfreq
if(mm > 1)
  line = fgetl(fid);
  strng = setstr(fgetl(fid));
  eval(['x = [',strng,'];']);
  FREQ = [FREQ x(1)];
  line = fgetl(fid);
end;
  for nn = 1:nphone
    strng = setstr(fgetl(fid));
    eval(['x = [',strng,'];']);
    depth(nn) = x;
  end;
  for nn = 1:nphone
    line = fgetl(fid);
    line = strrep(line,'(',' ');
    line = strrep(line,')',' ');
    line = strrep(line,',',' ');
    eval(['x = [',line,'];']);
    prr(nn,mm) = x(2)+i*x(3);
  end;
end;
saga_depth = depth;
saga_pr = prr;


