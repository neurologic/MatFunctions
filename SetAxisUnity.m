function SetAxisUnity(hfig)

figure(hfig)
axis tight
axis square
ylims = get(gca,'YLim');
xlims = get(gca,'XLim');

min_figure = min([xlims,ylims]);
max_figure = max([xlims,ylims]);

set(gca,'XLim',[min_figure,max_figure],'YLim',[min_figure,max_figure])
line([min_figure,max_figure],[min_figure,max_figure])