function hf = vis_colorbar(cmap)

hf = figure('Units','normalized'); 
colormap(cmap)
hCB = colorbar('north');
set(gca,'Visible',false);
hCB.Position = [0.15 0.3 0.74 0.4];
hCB.Ticks = [];
hCB.Color = [1 1 1 1];
hf.Position(4) = 0.1000;
hf.Units = 'pixels';
set(gcf,'color','w');

end