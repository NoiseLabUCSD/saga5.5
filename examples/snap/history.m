load -ascii sspmisab.mat 
 x=sspmisab(:,4:size(sspmisab,2)); 
 subplot(2,1,1)
plot (x(1:10,:)')
axis([1 4 0 51])
 subplot(2,1,2)  
 plot (x(1:100,:)') 
 axis([1 4 0 51])
 xlabel('Parameter')

 
 [mtms1,mps1,mtime1,mrng1]= ...
plotfield('./multstat.trf','timer',100,56.9999,'k-');

plot(mtime1*1000,mtms1(1,:),'r-','linewidth',1.5);
