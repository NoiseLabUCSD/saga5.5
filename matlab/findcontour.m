function ZHPD = findcontour1(zcontour,HPD);
totsum = sum(zcontour(:));
maxc = max(zcontour(:))*0.95;
del  = abs(max(zcontour(:)))./50;
zz   = fliplr([0:del:maxc]);
kk = 1;
%zcontour = zcontr;
for jj = 1:length(zz)
   [I] = find(zcontour(:) > zz(jj));
   frac(jj) = sum(zcontour(I))./totsum;
end
for jj = 1:length(HPD)
   hpd = HPD(jj);
   iz = max(find(frac.*100 < hpd));
   ZHPD(jj) = interp1([frac(iz) frac(iz+1)]*100,[zz(iz) zz(iz+1)],hpd,'spline');
end
