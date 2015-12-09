function runsaga(input,forw)
%input='horiz'
inputfil=input;
%forw='snap cov'

unix([' /bin/rm ',  input '.ext']);
[s,w]=unix([' saga ',input,' ',forw,'>list  &'])

%
% wait for saga to write first population
[a,b]=unix([' wc ' input '.ext']);
bsize=sscanf(b(6:10),'%i',1);
ti = clock; 

while (isempty(bsize)==1 | bsize==0);
  pause(10);
  [a,b]=unix([' wc ' input '.ext']);
  bsize=sscanf(b(6:10),'%i',1);
end  
timewait=etime(clock,ti);
%
% unix(ps -elf saga) 
while ~(unix ('top 4 |grep saga >/tmp/dummy')); 
  unix(['post ' input ' ' forw])
  %  plot_tc1v6lay
  plotsaga
  sagascat
  pause(timewait)
end

