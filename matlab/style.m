function style(taillefont,lw)
if nargin==0
   taillefont=11;
   lw=1.5;
 elseif nargin==1
   lw=2;
 end 
fig=gcf;
axes=get(fig,'currentaxes');
set(axes,'fontsize',taillefont);
set(axes,'LineWidth',1.5);
xla=get(axes,'xlabel');
yla=get(axes,'ylabel');
zla=get(axes,'zlabel');
ligne=get(axes,'children');
for I=1:length(ligne)
  if strcmp(get(ligne(I),'type'),'line')==1
    set(ligne(I),'linewidth',lw);
  end;
  if strcmp(get(ligne,'type'),'fontsize')==1
    set(ligne,'fontsize',taillefont);
  end;
end

set(xla,'fontsize',taillefont);
set(yla,'fontsize',taillefont);
set(zla,'fontsize',taillefont);

if nargin == 0
   fig=gcf;
end;
set(fig,'color',[1 1 1]);
