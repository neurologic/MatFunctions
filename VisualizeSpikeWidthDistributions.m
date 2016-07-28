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
isNarrow_TP = [];
isNarrow_HH = [];
width_ofmean_half_height = [];
width_ofmean_trough_peak = [];

allwavs = [];
for id = 1:size(sites,2)
% for id = 1:size(sites,2)
    load([sites{id} '/waveform_widths.mat'])
    thesesites = zeros(1,size(struct.cluster,2));
    thesesites(:) = id;
    siteID = [siteID, thesesites];
    clusterID = [clusterID, struct.cluster];
    isNarrow_TP = [isNarrow_TP, struct.narrow_trough_peak];
    isNarrow_HH = [isNarrow_HH, struct.narrow_half_height];
    isNarrow = [isNarrow, (struct.narrow_trough_peak .* struct.narrow_half_height)];
    width_ofmean_half_height = [width_ofmean_half_height struct.width_ofmean_half_height*1000];
    width_ofmean_trough_peak = [width_ofmean_trough_peak struct.width_ofmean_trough_peak*1000];
    
    load([sites{id} '/waveform_means.mat'])
    allwavs = [allwavs,arr];
end

load('../../MetaMat/metatoes.mat','allsortclass')

edges = [0:0.02:1];
n_halfheight = histc(width_ofmean_half_height,edges);
n_troughpeak = histc(width_ofmean_trough_peak,edges);

figure;hold on
scatter(width_ofmean_half_height(find(isNarrow==1)),width_ofmean_trough_peak(find(isNarrow==1)),100,'r','fill')
scatter(width_ofmean_half_height(find(isNarrow==0)),width_ofmean_trough_peak(find(isNarrow==0)),100,'b','fill')
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

datamat = [width_ofmean_half_height',width_ofmean_trough_peak'];

wavsmat = [];
for iwav = 1:size(allwavs,2)
    wavsmat(iwav,:) = allwavs{iwav}(1:120);
end
narrow_good = intersect(find(allsortclass==2),find(isNarrow==1));
wide_good = intersect(find(allsortclass==2),find(isNarrow==0));

figure;
hold on
plot(wavsmat(wide_good,:)','color','b')
plot(wavsmat(narrow_good,:)','color','r')


figure;
hold on
xtime = [1:size(wavsmat,2)]/20000*1000;
line(xtime,nanmean(wavsmat(narrow_good,:))','color','r','LineWidth',4)
line(xtime,nanmean(wavsmat(wide_good,:))','color','b','LineWidth',4)
set(gca,'YTick',[])
xlabel('msec')


% idx = kmeans(datamat,2);
% 
% figure;hold on
% scatter(width_ofmean_half_height(idx==1),width_ofmean_trough_peak(idx==1),100,'b','fill')
% scatter(width_ofmean_half_height(idx==2),width_ofmean_trough_peak(idx==2),100,'r','fill')
% xlabel('half-height width')
% ylabel('trough-to-peak width')

