hm = HeatMap(dist)
ax = hm.plot;
colormap('jet');
colorbar('Peer', ax)
caxis(ax, [0, max(dist(:))]);
set(ax, 'Xtick', 1:length(dist));
set(ax, 'Ytick', 1:length(dist));