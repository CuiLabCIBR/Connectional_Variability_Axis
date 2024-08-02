function plot_net_mat(net_mat, net_label, net_order)

t = [];
start = 1;
lines = 1;
mask = net_label;

for i = 1:length(net_order)
    add = find(mask==net_order(i));
    t = [t;add];
    start = start + length(add);
    lines(i+1) = start;
end

idx_begin = lines(1:end-1);
idx_end = lines(2:end)-1;
for i = 1:length(lines)-1
    for j = 1:length(lines)-1
        data_mat(idx_begin(i):idx_end(i),idx_begin(j):idx_end(j)) = net_mat(net_order(i),net_order(j));
    end
end

figure;
imagesc(data_mat);

hold on;
line_end = lines(end)-0.5;
for j = 2:length(lines) % draw lines dividing network
    line([0.5,line_end],[lines(j)-0.5,lines(j)-0.5],'Color','black','Linewidth',0.5);
    line([lines(j)-0.5,lines(j)-0.5],[0.5,line_end],'Color','black','Linewidth',0.5);
end
line([0.5,line_end],[0.5,0.5],'Color','black','Linewidth',0.5);

load('cmap.mat','cmap')
colormap(cmap);

set(gca,'xtick',[],'xticklabel',[])
set(gca,'ytick',[],'yticklabel',[])
box off
axis square;

set(gcf, 'units', 'inches', 'position', [0, 0, 5, 5], 'PaperUnits', 'inches', 'PaperSize', [5, 5])