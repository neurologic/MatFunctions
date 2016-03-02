%%%%%%%%%%%%%%%%%%
% this script is run from the MetaMat folder for each bird on Ice, which
% contains metatoes.mat and IntracellularData.mat

% RETURN:
% allSpikes
% IntracellularData
% potentially_not_auditory

%%
addpath('/Users/kperks/mnt/cube/Ice/kperks/MatFunctions/');
%% set default plot params
set(0,'DefaultAxesFontSize',20)
set(0,'DefaultTextFontSize',20)
set(0,'DefaultTextFontName','Helvetica')
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesColorOrder',[0 0 0;1 0 1;0 1 1; 0 1 0])

%% load up intracellular data
IntracellularData = load('IntracellularData.mat');
IntracellularData = IntracellularData.IntracellularData;

nsites_I = size(IntracellularData,2);

%% load up extracellular data
metatoes = load('metatoes.mat');

allsortclass = metatoes.allsortclass;
allchanz = metatoes.allchanz;
allsiteID = metatoes.allsiteID;

metatoes = metatoes.metatoes;

siteID_E = unique(allsiteID);
nsites_E = size(siteID_E,2);

% get spikes mat filtered with synaptic
trial_edges = [-2,8.5];
exptdt = IntracellularData{1}.dt;

type = 'synaptic';
trialnsamps = round(diff(trial_edges)*(1/exptdt));
%  trialnsamps = 10000;
xtime_spikes = linspace(trial_edges(1),trial_edges(2),trialnsamps);
allSpikes = SpikeMat_filtered(metatoes,type,trial_edges,trialnsamps);
% allSpikes is a cell array of function handles to pass xtime_spikes to for
% individual trials to build sub-matrices to analyze
%% separate allSpikes into prestim and peristim to send to brad
fs = 31250;
stimfs = 40000;
stimdur = (metatoes{isite}.stims{istim}.stim_end_times(1)/fs)-(metatoes{isite}.stims{istim}.stim_start_times(1)/fs);
stimstart = min(find(xtime_spikes >0));
stimstop = min(find(xtime_spikes >stimdur));
allSpikes_prestim = allSpikes(:,:,:,1:stimstart-1);
allSpikes_peristim = allSpikes(:,:,:,stimstart:stimstop);

%% which clusters are auditory?
isauditory = nan(size(allSpikes,1),size(allSpikes,2));
for isite = 1:size(metatoes,1)
    parfor istim = 1:size(metatoes{1}.stims,1)
        
        fs = 31250;
        stimfs = 40000;
        stimdur = (metatoes{isite}.stims{istim}.stim_end_times(1)/fs)-(metatoes{isite}.stims{istim}.stim_start_times(1)/fs);

        stimstart = min(find(xtime_spikes >0));
        stimstop = min(find(xtime_spikes >stimdur));
        
        
        ntrials = size(allSpikes,3);
        this_cluster_stim_response = zeros(ntrials,trialnsamps);
        for itrial = 1:ntrials
            this_cluster_stim_fun = allSpikes{isite,istim,itrial};
            if ~isempty(this_cluster_stim_fun);
                this_cluster_stim_response(itrial,:) = this_cluster_stim_fun(xtime_spikes);
            end
        end
        trial_avg = nanmean(this_cluster_stim_response,1);
        
        prestim = trial_avg(1:stimstart-1);
        peristim = trial_avg(stimstart:stimstop);
        
        [p,h,stats] = ranksum(prestim,peristim);
        % a "1" rejects the null that the prestim and peristim are from
        % continuous distributions with equal medians
        isauditory(isite,istim) = h;
    end
    isite
end

percent_driving_stims = 100*(sum(isauditory,2) ./ size(metatoes{1}.stims,1));
potentially_not_auditory = find(percent_driving_stims < 25)

%% relationship between mean signal correlation and primary recording site location
peristim_c = nan(size(metatoes,1),size(metatoes,1),size(metatoes{1}.stims,1));
prestim_c = nan(size(metatoes,1),size(metatoes,1),size(metatoes{1}.stims,1));
cluster_distance = [];

fs = 31250;
stimfs = 40000;
stimdur = (metatoes{1}.stims{1}.stim_end_times(1)/fs)-(metatoes{1}.stims{1}.stim_start_times(1)/fs);

stimstart = min(find(xtime_spikes >0));
stimstop = min(find(xtime_spikes >stimdur));

for istim = 1:size(metatoes{1}.stims,1)
    for isite1 = 1:size(metatoes,1)
        
        this_cluster_stim = allSpikes(isite1,istim,:,:);
        trial_avg_1 = squeeze(nanmean(this_cluster_stim,3));
        prestim_1 = trial_avg_1(1:stimstart-1);
        peristim_1 = trial_avg_1(stimstart:stimstop);
        
        for isite2 = 1:size(metatoes,1)
            
            this_cluster_stim = allSpikes(isite2,istim,:,:);
            trial_avg_2 = squeeze(nanmean(this_cluster_stim,3));
            prestim_2 = trial_avg_2(1:stimstart-1);
            peristim_2 = trial_avg_2(stimstart:stimstop);
            
            cluster_distance(isite1,isite2) = abs(allchanz(isite1) - allchanz(isite2));

            c = corrcoef(prestim_1,prestim_2);
            c = c(1,2);
            prestim_c(isite1,isite2,istim) = c;
            
            c = corrcoef(peristim_1,peristim_2);
            c = c(1,2);
            peristim_c(isite1,isite2,istim) = c;
        end
    end
end

mean_peristim_c = squeeze(nanmean(peristim_c,3));
mean_prestim_c = squeeze(nanmean(prestim_c,3));
a = ones(size(metatoes,1),size(metatoes,1));
inds = ~triu(a);
figure;
scatter(mean_peristim_c(inds),cluster_distance(inds))
xlabel('mean signal correlation (during stimulus)')
ylabel('primary site distance')
set(gca,'XLim',[-0.4,1])
text(0.4,1000,['n = ' num2str(size(metatoes,1)) ' clusters'])

figure;
scatter(mean_prestim_c(inds),cluster_distance(inds))
xlabel('mean signal correlation (during pre-stimulus silence)')
ylabel('primary site distance')
set(gca,'XLim',[-0.4,1])
text(0.4,1600,['n = ' num2str(size(metatoes,1)) ' clusters'])

%% spiking population stimulus response trial-averaged
prestim_dist_spiking = [];
peristim_dist_spiking = [];

n_prestim_spikes = [];
n_peristim_spikes = [];

siteref = 1;
for istim = 1:size(metatoes{1}.stims,1)
    fs = 31250;
    stimfs = 40000;
    stimdur = (metatoes{siteref}.stims{istim}.stim_end_times(1)/fs)-(metatoes{siteref}.stims{istim}.stim_start_times(1)/fs);
    
    stimstart = min(find(xtime_spikes >0));
    stimstop = min(find(xtime_spikes >stimdur));
    
    
    this_clusters_stims = allSpikes(:,istim,:,:);
    trial_avg = nanmean(this_clusters_stims,3);
    cluster_avg = squeeze(sum(trial_avg,1));
    cluster_avg = cluster_avg ./ max(cluster_avg);
    
    prestim_dist_spiking(istim,:) = cluster_avg(1:stimstart-1);
    peristim_dist_spiking(istim,:) = cluster_avg(stimstart:stimstop);
    

end

%% distribution spiking population stimulus response trial-averaged

edges = [-1:0.05:1];
n_prestim_spikes = histc(prestim_dist_spiking,edges,1)./size(prestim_dist_spiking,2);
n_peristim_spikes = histc(peristim_dist_spiking,edges,1)./size(prestim_dist_spiking,2);
figure;
hold on
stairs(edges, mean(n_peristim_spikes,2),'color','r','LineWidth',4)
stairs(edges, mean(n_prestim_spikes,2),'color','k','LineWidth',4)

yhist_peri_spikes = mean(n_peristim_spikes,2);
yhist_pre_spikes = mean(n_prestim_spikes,2);

%% CUBICm estimation on extracellular net population response
istim = 1;
S = n_peristim_spikes(istim,:);

tau = 10; %msec
% dt_spk = diff(trial_edges)/trialnsamps;
xihat = CuBICm(S,exptdt,{'exponential',1,tau/1000});
%% distribution intracellular stimulus response trial-averaged
prestim_dist_vm = [];
peristim_dist_vm = [];
edges = [-1:0.05:1];
   
    n_prestim_vm = [];
    n_peristim_vm = [];
figure;
hold on
for icell = 1:size(IntracellularData,2)
    thiscell = IntracellularData{icell};
 
    for istim = 1:size(thiscell.stim_data,2)
        
        normfactor = max([max(max(thiscell.prestim_data{istim})), max(max(thiscell.stim_data{istim}))]);
        prestim_dist_vm{istim} = mean(thiscell.prestim_data{istim})./normfactor;
        peristim_dist_vm{istim} = mean(thiscell.stim_data{istim})./normfactor;
        
        n_prestim_vm(icell,istim,:) = histc(prestim_dist_vm{istim},edges)./size(prestim_dist_vm{istim},2);
        n_peristim_vm(icell,istim,:) = histc(peristim_dist_vm{istim},edges)./size(prestim_dist_vm{istim},2);
    end
    
    stairs(edges, squeeze(mean(n_peristim_vm(icell,:,:),2)),'color','r','LineWidth',4)
    stairs(edges, squeeze(mean(n_prestim_vm(icell,:,:),2)),'color','k','LineWidth',4)
    
end

figure;
hold on
stairs(edges, squeeze(mean(mean(n_peristim_vm(icell,:,:),2),1)),'color','r','LineWidth',4)
stairs(edges, squeeze(mean(mean(n_prestim_vm(icell,:,:),2),1)),'color','k','LineWidth',4)

yhist_peri_vm = squeeze(mean(mean(n_peristim_vm(icell,:,:),2),1));
yhist_pre_vm = squeeze(mean(mean(n_prestim_vm(icell,:,:),2),1));

%% make log logs plot of distributions

log_edges = log(edges);
log_peri_spikes = log(yhist_peri_spikes);
log_peri_vm = log(yhist_peri_vm);

figure;
stairs(log_edges,log_peri_spikes,'color','m')
stairs(log_edges,log_peri_vm,'color','b')

%% stimulus dependence of correlations (able to split between noise and signal?)
%dependence of noise correlations on geometric mean between responses of neurons:

%%
edges = [min(allchanz):50:max(allchanz)];

%all clusters distribution
n = histc(allchanz,edges)
figure;
stairs(edges,n/size(allchanz,1),'LineWidth',4);
title('location distribution all clusters')
set(gca,'YLim',[0,0.15])
xlabel('depth in microns')

%good clusters distribution
n = histc(allchanz(find(allsortclass == 2)),edges);
figure;
stairs(edges,n/size(find(allsortclass == 2),1),'LineWidth',4)
title('location distribution good clusters')
set(gca,'YLim',[0,0.15])
xlabel('depth in microns')

%good clusters distribution
n = histc(allchanz(find(allsortclass == 1)),edges);
figure;
stairs(edges,n/size(find(allsortclass == 1),1),'LineWidth',4)
title('location distribution mua clusters')
set(gca,'YLim',[0,0.15])
xlabel('depth in microns')

%unsorted clusters distribution
n = histc(allchanz(find(allsortclass == 3)),edges);
figure;
stairs(edges,n/size(find(allsortclass == 1),1),'LineWidth',4)
title('location distribution mua clusters')
set(gca,'YLim',[0,0.15])
xlabel('depth in microns')


%
%
%% trial-averaged responses by depth (and population avg) for all stims their own figure
% normalize max cluster response for each cluster to 1
norm_allSpikes = [];
for icluster = 1:size(allSpikes,1)
    thiscluster = allSpikes(icluster,:,:,:);
    norm_scale = max(thiscluster(:));
    thiscluster = thiscluster./norm_scale;
    norm_allSpikes(icluster,:,:,:) = thiscluster;
end

maxplot_ht = 200;
scaled_Spikes = norm_allSpikes.*maxplot_ht;

include_inds = find(allsortclass);

% include_inds = find(allsortclass == 3);

for istim = 1:size(metatoes{1}.stims,1)
    figure;
    set(gcf,'Position',[680 83 449 1015])
    hold on
    
    xtime = linspace(-2,8.5,nsamps);
    for isite = 1:size(include_inds,1)
        
        dcoff = allchanz(include_inds(isite));
        
        meanpsth = mean(squeeze(scaled_Spikes(include_inds(isite),istim,:,:))) + dcoff;
        
        line(xtime,meanpsth')
        
    end
    this_clusters_stims = allSpikes(include_inds,istim,:,:);
    trial_avg = nanmean(this_clusters_stims,3);
    cluster_avg = sum(trial_avg,1);
    cluster_avg = cluster_avg ./ max(cluster_avg);
    line(xtime,((200*squeeze(cluster_avg))+max(allchanz)+100)','LineWidth',1)
    axis tight
    
    fs = 31250;
    stimfs = 40000;
    stimdur = (metatoes{1}.stims{1}.stim_end_times(1)/fs)-(metatoes{1}.stims{1}.stim_start_times(1)/fs);
    line([0,0],[min(allchanz),max(allchanz)])
    line([stimdur/6,stimdur/6],[min(allchanz),max(allchanz)])
    line([stimdur,stimdur],[min(allchanz),max(allchanz)])
end

%% for a single trial


%% population spiking magnitude distribution (trial-averaged psth for each stimulus)
% ALL CLUSTERS
include_inds = find(allsortclass);

sort_ind = 2;
% include_inds = find(allsortclass == sort_ind);

stimdur = (metatoes{1}.stims{1}.stim_end_times(1)/fs)-(metatoes{1}.stims{1}.stim_start_times(1)/fs);
xtime = linspace(-2,8.5,nsamps);
stim_onset_ind = crossing(xtime);
stim_offset_ind = crossing(xtime-stimdur);


meanpsth_across_clusters = [];
for istim = 1:16
    meanpsth_across_trials = squeeze(mean(allSpikes(include_inds,istim,:,stim_onset_ind:stim_offset_ind),3));
    meanpsth_across_clusters(istim,:) = mean(meanpsth_across_trials,1);
end
normfactor = max(max(meanpsth_across_clusters));
% meanpsth_across_clusters = meanpsth_across_clusters ./normfactor;
edges = [0:normfactor/100:normfactor];


distributions = [];
for istim = 1:16
    n = histc(meanpsth_across_clusters(istim,:),edges);
    distributions(istim,:) = n ./ size(meanpsth_across_clusters,2);
end
figure;
line(edges,distributions')

title(['only cluster' num2str(sort_ind)]);
%% for each site (so simultaneous data only) population spiking magnitude distribution (trial averaged and residuals on individual trials)

% use toedata3 as example first
this_toe = toedata4;

toedata_singleSite_psth = [];
spike_stimorder_SingleSite = [];

tau = 0.01;
nsamps =1000;

for isite = 1:size(this_toe,1)
    sortclass_singleSite(isite) = this_toe{isite}.sort_class;
    thissite = this_toe{isite};
    
    
    for istim = 1:size(this_toe{isite}.stims,1);
        thisstim = thissite.stims{istim};
        spike_stimorder_SingleSite{istim} = thisstim.name{:};
        for itrial = 1:thisstim.ntrials
            
            thisstim_toes = thisstim.toes;
            
            thistrial_toes = thisstim_toes{itrial};
            
            filteredspikes = zeros(1,nsamps);
            
            if ~isempty(thistrial_toes)
                fun = filtered_response_synaptic(thistrial_toes, tau);
                filteredspikes = fun(linspace(-2,8.5,nsamps));
            end
            
            toedata_singleSite_psth(isite,istim,itrial,:) = filteredspikes;
        end
    end
end

a = cell2mat(spike_stimorder_SingleSite');
spike_stimorder = a(:,1)';
%% toedata3 trial-avg distributions
include_inds = find(sortclass_singleSite);
% include_inds = find(allsortclass == 3);

meanpsth_across_clusters = [];
for istim = 1:16
    meanpsth_across_trials = squeeze(mean(toedata_singleSite_psth(include_inds,istim,:,stim_onset_ind:stim_offset_ind),3));
    meanpsth_across_clusters(istim,:) = mean(meanpsth_across_trials,1);
end
normfactor = max(max(meanpsth_across_clusters));
% meanpsth_across_clusters = meanpsth_across_clusters ./normfactor;
edges = [0:normfactor/100:normfactor];

stimdur = (metatoes{1}.stims{1}.stim_end_times(1)/fs)-(metatoes{1}.stims{1}.stim_start_times(1)/fs);
xtime = linspace(-2,8.5,nsamps);
stim_onset_ind = crossing(xtime);
stim_offset_ind = crossing(xtime-stimdur);

distributions = [];
for istim = 1:16
    n = histc(meanpsth_across_clusters(istim,:),edges);
    distributions(istim,:) = n ./ size(meanpsth_across_clusters,2);
end
figure;
line(edges,distributions')

title('Site3; all clusters; trial-avg psth; stimulus period')
%% toedata3 trial #10 distributions
include_inds = find(sortclass_singleSite);
% include_inds = find(allsortclass == 3);
%
% itrial = 20;
meanpsth_across_clusters = [];
for itrial = 1:20
    for istim = 1:16
        meanpsth_single_trials = squeeze(toedata_singleSite_psth(include_inds,istim,itrial,stim_onset_ind:stim_offset_ind));
        meanpsth_across_clusters(istim,itrial,:) = mean(meanpsth_single_trials,1);
    end
end
normfactor = max(max(max(meanpsth_across_clusters)));
% meanpsth_across_clusters = meanpsth_across_clusters ./normfactor;
edges = [0:normfactor/100:normfactor];

stimdur = (metatoes{1}.stims{1}.stim_end_times(1)/fs)-(metatoes{1}.stims{1}.stim_start_times(1)/fs);
xtime = linspace(-2,8.5,nsamps);
stim_onset_ind = crossing(xtime);
stim_offset_ind = crossing(xtime-stimdur);

distributions = [];
for istim = 1:16
    thisstim_alltrials = squeeze(meanpsth_across_clusters(istim,:,:));
    n = histc(thisstim_alltrials(:),edges);
    distributions(istim,:) = n ./ size(thisstim_alltrials(:),1);
end
figure;
line(edges,distributions')

title('Site3; all clusters; average single-trial dist')

%% for all stim, plot psth for all trials in black and trial-averaged psth in red
include_inds = find(sortclass_singleSite); % for all clusters in toedata3;
nstim = 16;
ncol = 2;
figure;
hold on
for istim = 1:nstim
    % istim = 16;
    
    subplot(ceil(nstim/ncol), ncol, istim)
    meanpsth_across_clusters = [];
    for itrial = 1:20
        %     for istim = 1:16
        meanpsth_across_trials = squeeze(toedata_singleSite_psth(include_inds,istim,itrial,:));
        meanpsth_across_clusters(istim,itrial,:) = mean(meanpsth_across_trials,1);
        line(xtime,mean(meanpsth_across_trials,1),'color','k');
        %     end
    end
    line(xtime,squeeze(mean(meanpsth_across_clusters,2)),'color','r')
    axis tight
    
    if istim ~= nstim
        set(gca,'XTickLabel',[])
    end
end

%% load('intracell_preprocessed.mat')

ICdata = intracell_preprocessed{3};
nstim = size(ICdata.stimnames,2);
ncol = 2;
figure;
hold on
for istim = 1:nstim
    thistrial_full = [ICdata.prestim_data{istim} ICdata.stim_data{istim} ICdata.poststim_data{istim}];
    subplot(ceil(nstim/ncol), ncol, istim)
    xtime = [1:size(thistrial_full,2)]*ICdata.dt;
    xtime_lin = linspace(-(ICdata.stimtimes(1)*ICdata.dt),(size(thistrial_full,2)*ICdata.dt)-(ICdata.stimtimes(1)*ICdata.dt),size(thistrial_full,2));
    line(xtime_lin,thistrial_full,'color','k');
    
    line(xtime_lin,(mean(thistrial_full,1)),'color','r')
    axis tight
    
    if istim ~= nstim
        set(gca,'XTickLabel',[])
    end
end


%% for all stim, plot psth for all trials in black and trial-averaged psth in red
include_inds = find(sortclass_singleSite); % for all clusters in toedata3;
intracellID = 4;
nstim = 8;
ncol = 2;

figure;
hold on
for istim = 1:nstim
   hs(istim) =  subplot(ceil(nstim/ncol), ncol, istim);
   hold on
end
for istim = 1:nstim
 
    ICdata = intracell_preprocessed{intracellID};
    stimname = ICdata.stimnames{istim};
    xtime_vm = linspace(-(ICdata.stimtimes(1)*ICdata.dt),(size(thistrial_full,2)*ICdata.dt)-(ICdata.stimtimes(1)*ICdata.dt),size(thistrial_full,2));
    thistrial_full = [ICdata.prestim_data{istim} ICdata.stim_data{istim} ICdata.poststim_data{istim}];
    trialavg_vm = mean(thistrial_full,1);
    trialavg_vm_norm = trialavg_vm ./ max(trialavg_vm);
    
    spike_stim_ind = regexp(spike_stimorder,stimname);
    meanpsth_across_clusters = [];
    for itrial = 1:20
        meanpsth_across_trials = squeeze(toedata_singleSite_psth(include_inds,spike_stim_ind,itrial,:));
        meanpsth_across_clusters(itrial,:) = mean(meanpsth_across_trials,1);
    end
    xtime_spike = linspace(-2,8.5,nsamps); %nsamps was nsamps used for linear filter
    trialavg_spike = mean(meanpsth_across_clusters,1);
    trialavg_spike_norm = trialavg_spike ./ max(trialavg_spike);
    
    hl1=line(xtime_vm,trialavg_vm_norm,'Color','k');
    hold on
    ax(1)=hs(istim);
    set(ax(1),'XColor','k','YColor','k');
    
    ax(2)=axes('Position',get(ax(1),'Position'),...
        'XAxisLocation','top',...
        'YAxisLocation','right',...
        'Color','none',...
        'XColor','r','YColor','k');
        
    hl2=line(xtime_spike,trialavg_spike_norm,'Color','r','Parent',ax(2));
    
    set(ax(1),'XLim',[0,6])
    ylabel(ax(1),stimname)
    set(ax(2),'XLim',[0,6.615])
    ylabel(ax(2),spike_stimorder(spike_stim_ind))
    set(ax,'YLim',[-0.5,1])

    
    set(ax,'YTickLabel',[])
    if istim ~= nstim
        set(ax,'XTickLabel',[])
    end
    
end
%% for all metatoes do this same plot...
%for all stim, plot psth for all trials in black and trial-averaged psth in red
include_inds = find(allsortclass); % for all clusters in toedata3;
nstim = 16;
ncol = 2;
figure;
hold on
for istim = 1:nstim
    % istim = 16;
    
    subplot(ceil(nstim/ncol), ncol, istim)
    meanpsth_across_clusters = [];
    for itrial = 1:20
        %     for istim = 1:16
        meanpsth_across_trials = squeeze(allSpikes(include_inds,istim,itrial,:));
        meanpsth_across_clusters(itrial,:) = mean(meanpsth_across_trials,1);
        %         line(xtime,mean(meanpsth_across_trials,1),'color','k');
        %%%%%%%%%%%%% since this is across sites, plotting each trial seems uninformative)
        %%%%%%%%%%%% only comparing trial-averaged mean seems informative%     end
    end
    line(xtime,mean(meanpsth_across_clusters,1),'color','r')
    axis tight
    
    if istim ~= nstim
        set(gca,'XTickLabel',[])
    end
end

