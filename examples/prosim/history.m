
[mtms1,mps1,mtime1,mrng1]= ...
plotfield('./prosim_out.dat','timer',100,56.9999,'k-');

plot(mtime1*1000,mtms1(1,:),'r-','linewidth',1.5);
