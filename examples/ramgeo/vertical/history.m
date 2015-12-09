tc='3'
casedir=[ '/heart5/gerstoft/itworkshop/tc' tc '/']
!/bin/rm pres.in
fr=[25:200 205:5:500];
for freq=fr;
if (freq<100)
  [pr par]  =read_it([casedir 'ivwkt' tc '_v_00' num2str(freq) '.cpr']);
else
  [pr par]  =read_it([casedir 'ivwkt' tc '_v_0' num2str(freq) '.cpr']);
end  
depth=par(5)+par(6)*[0:(par(7)-1)];
write_pvec(pr(:,5),freq,depth)
end
