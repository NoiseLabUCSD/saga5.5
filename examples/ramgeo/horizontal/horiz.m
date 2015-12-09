tc='0'
casedir=[ '/net/heart/heart5/gerstoft/itworkshop/tc' tc '/']
!/bin/rm pres.in
fr=[25:200 205:5:500];
%fr=[100 200 300 400 500];
%fr=200
ifreq=0;
prbot=[];
for freq=fr;
ifreq=ifreq+1;
  if (freq<100)
  [pr par]  =read_it([casedir 'ivwkt' tc '_h_00' num2str(freq) '.cpr']);
else
  [pr par]  =read_it([casedir 'ivwkt' tc '_h_0' num2str(freq) '.cpr']);
end  
depth=par(5)+par(6)*[0:(par(7)-1)];
range=par(2)+par(3)*[0:(par(4)-1)];
%write_pvec(pr(:,5),freq,depth)
istart=200*1;  %1000 m 
iend=istart+60;    %1000+300 m
write_pvec(pr(2,istart:iend).',freq,range(istart:iend))
prbot(ifreq,:)=pr(2,:);
end

iloop=1
istart=200*iloop;
iend=istart+60;
!rm pres.in
write_dformat(prbot(:,istart:(iend)),'horizontal arrray 1000-1300')
