%VisualizeSpikeWidthDistributions
%%%%%%%%
%run from bird's phy directory so that don't have to define path
%%%%%%%%

% path_to_phy = 'Z:\Ice\kperks\B1112\klusta\phy041316\';
% path_to_phy = '..\klusta\phy041316\';
sites = [];
d = dir();
for id = 1:size(d,1)-2
    sites{id} = d(id+2).name;
end
clusterID = [];
siteID = [];
isNarrow = [];
widths_ofmean_half_height = [];
widths_ofmean_trough_peak = [];

for id = [1,2,3]
% for id = 1:size(sites,2)
    load([sites{id} '/waveform_widths.mat'])
    thesesites = zeros(1,size(struct.cluster,2));
    thesesites(:) = id;
    siteID = [siteID, thesesites];
    clusterID = [clusterID, struct.cluster];
    isNarrow = [isNarrow, struct.narrow_trough_peak];
    widths_ofmean_half_height = [widths_ofmean_half_height struct.widths_ofmean_half_height*1000];
    widths_ofmean_trough_peak = [widths_ofmean_trough_peak struct.widths_ofmean_trough_peak*1000];
end

edges = [0:0.02:1];
n_halfheight = histc(widths_ofmean_half_height,edges);
n_troughpeak = histc(widths_ofmean_trough_peak,edges);

figure;
scatter(widths_ofmean_half_height,widths_ofmean_trough_peak)
xlabel('half-height width')
ylabel('trough-to-peak width')

figure;
stairs(edges,n_halfheight)
xlabel('edges')
ylabel('half-height width')

figure;
stairs(edges,n_troughpeak)
xlabel('edges')
ylabel('trough-to-peak width')