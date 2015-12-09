  imagesc(x,y,rot90(c_val(:,:))), colormap(1-gray), ...
          shading('flat'),h=colorbar;  %caxis([-7 0]),
 set(gca,'YDir','normal')
xlabel('Range (km)')
ylabel('Height (m)')
title('Fig 5')
