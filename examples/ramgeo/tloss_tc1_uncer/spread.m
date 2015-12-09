function [sp ]=spread(no,xval)
% returns lower, upper and median fractile
% 
  no1=no/sum(no);
  perc=0.05;
  ss1=0; ii1=0; ii2=0; ii3=0;
 nsamp=length(xval);
 for ii=1:nsamp
    ss1=ss1+no1(ii);
    if (ss1>perc &ii1==0); ii1=ii; lsum=ss1; end
    if (ss1>1-perc & ii2==0); ii2=ii; usum=ss1; end
    if (ss1>0.5 & ii3==0); ii3=ii; msum=ss1; end
  end
%     spreadp(ir,id)=xval(ii2)-xval(ii1);   medp(ir,id)=xval(ii3);
if (ii1>1)
  sp(1)=xval(ii1)-(xval(ii1)-xval(ii1-1))/no1(ii1)*(lsum-perc); 
else
  sp(1)=xval(ii1);
end

%keyboard
if (ii2~=1 )
  sp(2)=xval(ii2)-(xval(ii2)-xval(ii2-1))/no1(ii2)*(usum-(1-perc)); 
else
  sp(2)=xval(ii2);
end

if (ii3>1)
  sp(3)=xval(ii3)-(xval(ii3)-xval(ii3-1))/no1(ii3)*(msum-0.5); 
else
  sp(3)=xval(ii3);
end



%keyboard
%end