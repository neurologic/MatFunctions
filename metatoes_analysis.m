%%%%%%%%%%%%%%%%%%
% this script is run from the MetaMat folder for each bird on Ice, which
% contains metatoes.mat and IntracellularData.mat

% RETURN:
% allSpikes
% IntracellularData
% potentially_not_auditory

kwik2mat_KP

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
metatoesVars = metatoes;
vars = fieldnames(metatoes);
for ifield = 1:size(vars,1)
    s = [vars{ifield} '= metatoesVars.' vars{ifield} ';'];
    eval(s)
end

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

fs = 31250;
stimfs = 40000;
stimdur = (metatoes{1}.stims{1}.stim_end_times(1)/fs)-(metatoes{1}.stims{1}.stim_start_times(1)/fs);
stimstart = min(find(xtime_spikes >0));
stimstop = min(find(xtime_spikes >stimdur));

clear metatoes
% 
% allSpikes_prestim = allSpikes(:,:,:,1:stimstart-1);
% allSpikes_peristim = allSpikes(:,:,:,stimstart:stimstop);

%% on lintu 
%which clusters are auditory?
% isauditory = nan(size(allSpikes,1),size(allSpikes,2));
% for icluster = 1:size(metatoes,1)
%     parfor istim = 1:size(metatoes{1}.stims,1)
%         
%         fs = 31250;
%         stimfs = 40000;
%         
%         ntrials = size(allSpikes,3);
%         this_cluster_stim_response = zeros(ntrials,trialnsamps);
%         for itrial = 1:ntrials
%             this_cluster_stim_fun = allSpikes{icluster,istim,itrial};
%             if ~isempty(this_cluster_stim_fun);
%                 this_cluster_stim_response(itrial,:) = this_cluster_stim_fun(xtime_spikes);
%             end
%         end
%         trial_avg = nanmean(this_cluster_stim_response,1);
%         
%         prestim = trial_avg(1:stimstart-1);
%         peristim = trial_avg(stimstart:stimstop);
%         
%         [p,h,stats] = ranksum(prestim,peristim);
%         % a "1" rejects the null that the prestim and peristim are from
%         % continuous distributions with equal medians
%         isauditory(icluster,istim) = h;
%     end
%     icluster
% end
%% which clusters are auditory?
isauditory = load('WorkingData/isauditory.mat');
isauditory = isauditory.isauditory;

percent_driving_stims = 100*(sum(isauditory,2) ./ size(allSpikes,2));
potentially_not_auditory = find(percent_driving_stims < 25)

%% relationship between mean signal correlation and primary recording site location
peristim_c = nan(size(metatoes,1),size(metatoes,1),size(metatoes{1}.stims,1));
prestim_c = nan(size(metatoes,1),size(metatoes,1),size(metatoes{1}.stims,1));
cluster_distance = [];

fs = 31250;
stimfs = 40000;


for istim = 1:size(metatoes{1}.stims,1)
    for icluster1 = 1:size(metatoes,1)
        
        this_cluster_stim = allSpikes(icluster1,istim,:,:);
        trial_avg_1 = squeeze(nanmean(this_cluster_stim,3));
        prestim_1 = trial_avg_1(1:stimstart-1);
        peristim_1 = trial_avg_1(stimstart:stimstop);
        
        for icluster2 = 1:size(metatoes,1)
            
            this_cluster_stim = allSpikes(icluster2,istim,:,:);
            trial_avg_2 = squeeze(nanmean(this_cluster_stim,3));
            prestim_2 = trial_avg_2(1:stimstart-1);
            peristim_2 = trial_avg_2(stimstart:stimstop);
            
            cluster_distance(icluster1,icluster2) = abs(allchanxz(icluster1) - allchanxz(icluster2));

            c = corrcoef(prestim_1,prestim_2);
            c = c(1,2);
            prestim_c(icluster1,icluster2,istim) = c;
            
            c = corrcoef(peristim_1,peristim_2);
            c = c(1,2);
            peristim_c(icluster1,icluster2,istim) = c;
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
ylabel('primary cluster distance')
set(gca,'XLim',[-0.4,1])
text(0.4,1000,['n = ' num2str(size(metatoes,1)) ' clusters'])

figure;
scatter(mean_prestim_c(inds),cluster_distance(inds))
xlabel('mean signal correlation (during pre-stimulus silence)')
ylabel('primary cluster distance')
set(gca,'XLim',[-0.4,1])
text(0.4,1600,['n = ' num2str(size(metatoes,1)) ' clusters'])


%% distribution spiking population stimulus response trial-averaged
% from processing done on lintu
this_stim_response = load('WorkingData/this_stim_response.mat');
this_stim_response = this_stim_response.this_stim_response;

prestim_dist_spiking = [];
peristim_dist_spiking = [];
prestim_skew = [];
prestim_var = [];
peristim_skew = [];
peristim_var = [];

for istim = 1:size(this_stim_response,1);
    pop_response = this_stim_response(istim,:);
    pop_response = pop_response - min(pop_response);
    pop_response = pop_response./(max(pop_response));
    
    prestim_dist_spiking(istim,:) = pop_response(1:stimstart-1);
    prestim_skew(istim,:) = skewness(pop_response(1:stimstart-1));
    prestim_var(istim,:) = var(pop_response(1:stimstart-1));
    prestim_mean(istim,:) = mean(pop_response(1:stimstart-1));
    
    peristim_dist_spiking(istim,:) = pop_response(stimstart:stimstop);
    peristim_skew(istim,:) = skewness(pop_response(stimstart:stimstop));
    peristim_var(istim,:) = var(pop_response(stimstart:stimstop));
    peristim_mean(istim,:) = mean(pop_response(stimstart:stimstop));
end

edges = [-1:0.05:1];
n_prestim_spikes = histc(prestim_dist_spiking,edges,1)./size(prestim_dist_spiking,2);
n_peristim_spikes = histc(peristim_dist_spiking,edges,1)./size(prestim_dist_spiking,2);
figure;
hold on
stairs(edges, mean(n_peristim_spikes,2),'color','r','LineWidth',4)
stairs(edges, mean(n_prestim_spikes,2),'color','k','LineWidth',4)

yhist_peri_spikes = mean(n_peristim_spikes,2);
yhist_pre_spikes = mean(n_prestim_spikes,2);

hfig = figure;
scatter(prestim_skew, peristim_skew, 100,'k','fill')
SetAxisUnity(hfig)
xlabel('pre-stimulus skewness spiking dist')
ylabel('peri-stim skewness spiking dist')

hfig = figure;
scatter(prestim_var, peristim_var, 100,'k','fill')
SetAxisUnity(hfig)
xlabel('pre-stimulus variance spiking dist')
ylabel('peri-stim variance spiking dist')

hfig = figure;
scatter(prestim_mean, peristim_mean, 100,'k','fill')
SetAxisUnity(hfig)
xlabel('pre-stimulus mean spiking dist')
ylabel('peri-stim mean spiking dist')

clear this_stim_response
%% distribution spiking population stimulus response single trials
% from processing done on lintu
this_trial_response = load('WorkingData/this_trial_response.mat');
this_trial_response = this_trial_response.this_trial_response;

prestim_dist_spiking = [];
peristim_dist_spiking = [];
prestim_skew = [];
prestim_var = [];
peristim_skew = [];
peristim_var = [];

for istim = 1:size(this_trial_response,1);
    for itrial = 1:size(this_trial_response,2);
        pop_response = squeeze(this_trial_response(istim,itrial,:));

        pop_response = pop_response - min(pop_response);
        pop_response = pop_response./(max(pop_response));
        
        prestim_response = pop_response(:,1:stimstart-1);
        prestim_dist_spiking(istim,:) = prestim_response(:);
        prestim_skew(istim,:) = skewness(prestim_response(:));
        prestim_var(istim,:) = var(prestim_response(:));
        prestim_mean(istim,:) = mean(prestim_response(:));
        
        peristim_response = pop_response(:,stimstart:stimstop);
        peristim_dist_spiking(istim,:) = peristim_response(:);
        peristim_skew(istim,:) = skewness(peristim_response(:));
        peristim_var(istim,:) = var(peristim_response(:));
        peristim_mean(istim,:) = mean(peristim_response(:));
    end
end

edges = [-1:0.05:1];
n_prestim_spikes = histc(prestim_dist_spiking,edges,1)./size(prestim_dist_spiking,2);
n_peristim_spikes = histc(peristim_dist_spiking,edges,1)./size(prestim_dist_spiking,2);
figure;
hold on
stairs(edges, mean(n_peristim_spikes,2),'color','r','LineWidth',4)
stairs(edges, mean(n_prestim_spikes,2),'color','k','LineWidth',4)

hfig = figure;
scatter(prestim_skew, peristim_skew, 100,'k','fill')
SetAxisUnity(hfig)
xlabel('pre-stimulus skewness spiking dist')
ylabel('peri-stim skewness spiking dist')

hfig = figure;
scatter(prestim_var, peristim_var, 100,'k','fill')
SetAxisUnity(hfig)
xlabel('pre-stimulus variance spiking dist')
ylabel('peri-stim variance spiking dist')

hfig = figure;
scatter(prestim_mean, peristim_mean, 100,'k','fill')
SetAxisUnity(hfig)
xlabel('pre-stimulus mean spiking dist')
ylabel('peri-stim mean spiking dist')

clear this_trial_response

%% spiking population stimulus response per trial



for isite = 1:nsite_E
    clusterinds = find(allsiteID == siteID_E(isite));

    this_trial_response = zeros(nstim,ntrials,trialnsamps);
    
    for istim = 1:nstim
        fs = 31250;
        stimfs = 40000;
        stimdur = (metatoes{clusterref}.stims{istim}.stim_end_times(1)/fs)-(metatoes{clusterref}.stims{istim}.stim_start_times(1)/fs);
        
        for itrial = 1:ntrials
            
            this_cluster_response = zeros(size(clusterinds,2),trialnsamps);
            for icluster = 1:size(clusterinds,2)
                this_cluster_stim_fun = allSpikes{clusterinds(icluster),istim,itrial};
                if ~isempty(this_cluster_stim_fun);
                    this_cluster_response(icluster,:) = this_cluster_stim_fun(xtime_spikes);
                end
            end
            trial_avg = sum(this_cluster_response,1);
            this_trial_response(istim,itrial,:) = trial_avg;
        end
        
    end
    
    save(['WorkingData/this_trial_response_site' num2str(isite) '.mat'],'this_trial_response')
end
%% CUBICm estimation on extracellular net population response
istim = 1;
S = peristim_dist_spiking(istim,:);

tau = 10; %msec
% dt_spk = diff(trial_edges)/trialnsamps;
xihat = CuBICm(S,exptdt,{'exponential',1,tau/1000});


% %% spiking population stimulus response per trial
% prestim_dist_spiking = [];
% peristim_dist_spiking = [];
% 
% n_prestim_spikes = [];
% n_peristim_spikes = [];
% 
% nstim = size(allSpikes,2);
% clusterref = 1;
% for istim = 1:nstim
%     fs = 31250;
%     stimfs = 40000;
%     stimdur = (metatoes{clusterref}.stims{istim}.stim_end_times(1)/fs)-(metatoes{clusterref}.stims{istim}.stim_start_times(1)/fs);
%     
%     
%     ntrials = size(allSpikes,3);
%     this_trial_response = zeros(ntrials,trialnsamps);
%     for itrial = 1:ntrials
%         
%         ncluster = size(allSpikes,1);
%         this_cluster_response = zeros(ncluster,trialnsamps);
%         for icluster = 1:ncluster
%             this_cluster_stim_fun = allSpikes{icluster,istim,itrial};
%             if ~isempty(this_cluster_stim_fun);
%                 this_cluster_response(icluster,:) = this_cluster_stim_fun(xtime_spikes);
%             end
%         end
%         trial_avg = sum(this_cluster_response,1);
%         this_trial_response(istim,trial,:) = trial_avg;
%     end
%    
% end
% 
% save('WorkingData/this_trial_response.mat','this_trial_response')
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
edges = [min(allchanxz):50:max(allchanxz)];

%all clusters distribution
n = histc(allchanxz,edges)
figure;
stairs(edges,n/size(allchanxz,1),'LineWidth',4);
title('location distribution all clusters')
set(gca,'YLim',[0,0.15])
xlabel('depth in microns')

%good clusters distribution
n = histc(allchanxz(find(allsortclass == 2)),edges);
figure;
stairs(edges,n/size(find(allsortclass == 2),1),'LineWidth',4)
title('location distribution good clusters')
set(gca,'YLim',[0,0.15])
xlabel('depth in microns')

%good clusters distribution
n = histc(allchanxz(find(allsortclass == 1)),edges);
figure;
stairs(edges,n/size(find(allsortclass == 1),1),'LineWidth',4)
title('location distribution mua clusters')
set(gca,'YLim',[0,0.15])
xlabel('depth in microns')

%unsorted clusters distribution
n = histc(allchanxz(find(allsortclass == 3)),edges);
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

% include_inds = find(allsortclass);

include_inds = find(allsortclass == 2);

ntrials = size(allSpikes_SynFilter,3);
nsamps = size(xtime_spikes,2);
for istim = 8%:size(metatoes{1}.stims,1)
    figure;
    set(gcf,'Position',[680 83 449 1015],'Visible','off')
    hold on
    
%     xtime = linspace(-2,8.5,nsamps);
    for icluster = 1:size(include_inds,1)
        
        dcoff = allchanxz(include_inds(icluster));
        
        ClusterFun = squeeze(allSpikes_SynFilter(include_inds(icluster),istim,:));
        
        ClusterFiltR = zeros(ntrials,nsamps);
        for itrial = 1:ntrials
            trialfun = ClusterFun{itrial};
            if ~isempty(trialfun)
                ClusterFiltR(itrial,:) = trialfun(xtime_spikes);
            end
        end
        meanpsth = mean(ClusterFiltR,1);
        meanpsth = ((meanpsth./max(meanpsth))*25) + dcoff;
        
        line(xtime_spikes,meanpsth')
        
    end
    set(gcf,'Visible','on')
    
end
    this_clusters_stims = allSpikes(include_inds,istim,:,:);
    trial_avg = nanmean(this_clusters_stims,3);
    cluster_avg = sum(trial_avg,1);
    cluster_avg = cluster_avg ./ max(cluster_avg);
    line(xtime,((200*squeeze(cluster_avg))+max(allchanxz)+100)','LineWidth',1)
    axis tight
    
    fs = 31250;
    stimfs = 40000;
    stimdur = (metatoes{1}.stims{1}.stim_end_times(1)/fs)-(metatoes{1}.stims{1}.stim_start_times(1)/fs);
    line([0,0],[min(allchanxz),max(allchanxz)])
    line([stimdur/6,stimdur/6],[min(allchanxz),max(allchanxz)])
    line([stimdur,stimdur],[min(allchanxz),max(allchanxz)])
end
%%
do_trialshuffle = 0;

%defaults inds
sortclassinds = [1:size(allsortclass,2)];
siteinds = [1:size(allsiteID,2)];

sortclassinds = find(allsortclass == 2);
siteinds = find(allsiteID == 3);
include_inds = intersect(sortclassinds,siteinds);

nclusters = size(include_inds,1);
nstims = size(allSpikes_SynFilter,2);
ntrials = size(allSpikes_SynFilter,3);
nsamps = size(xtime_spikes,2);

ResponseSumFiltR = zeros(nstims,ntrials,nsamps);

for istim = 1:size(metatoes{1}.stims,1)
    
    TrialFiltR = zeros(nclusters,nsamps);
    
%     hfig = figure;
%     hold on
       
    
    for icluster = 1:size(include_inds,1)
        trialshuffle(icluster,:) = randperm(ntrials);
    end
    
    FiltFun = squeeze(allSpikes_SynFilter(include_inds,istim,:));
    for itrial = 1:ntrials
        
        for icluster = 1:size(include_inds,1)

            trialind = itrial;
            if do_trialshuffle ==1
                trialind = trialshuffle(icluster,itrial);
            end
            
            clusterfun = FiltFun{icluster,trialind};
            if ~isempty(clusterfun)
                TrialFiltR(icluster,:) = clusterfun(xtime_spikes); %for this cluster on this trial
            end
        end
        meanpsth(itrial,:) = sum(TrialFiltR,1);
        
%         line(xtime_spikes,squeeze(meanpsth(itrial,:))')
        
    end
    
    ResponseSumFiltR(istim,:,:) = meanpsth;
    %     set(gcf,'Visible','on')
%     line(xtime_spikes,squeeze(mean(ResponseSumFiltR(istim,:,:),2)),'color','r','LineWidth',3);
%     set(gca,'YLim',[0,20]);

end


%%

ResidNorm = DataNorm(1,:,:) - repmat(mean(DataNorm(1,:,:),2),1,ntrials,1);
ResidShuffle = DataShuffle(1,:,:) - repmat(mean(DataShuffle(1,:,:),2),1,ntrials,1);

edges = [0:0.5:20];
bytrial = squeeze(ResidNorm(1,:,1:sigon));
skew_pre_norm = skewness(bytrial(:))
n_pre_norm = histc(bytrial(:),edges) / (sigon*ntrials);
bytrial = squeeze(ResidNorm(1,:,sigon:sigoff));
skew_stim_norm = skewness(bytrial(:))
n_stim_norm = histc(bytrial(:),edges) / ((sigoff-sigon)*ntrials);

edges = [0:0.5:20];
bytrial = squeeze(ResidShuffle(1,:,1:sigon));
skew_pre_shuffle = skewness(bytrial(:))
n_pre_shuffle = histc(bytrial(:),edges) / (sigon*ntrials);
bytrial = squeeze(ResidShuffle(1,:,sigon:sigoff));
skew_stim_shuffle = skewness(bytrial(:))
n_stim_shuffle = histc(bytrial(:),edges) / ((sigoff-sigon)*ntrials);

figure;hold on
stairs(edges,n_pre_norm,'color',[0,0,0],'LineWidth',5)
stairs(edges,n_pre_shuffle,'color',[0.8,0.8,0.8],'LineWidth',3)
stairs(edges,n_stim_norm,'color',[1,0,0],'LineWidth',5)
stairs(edges,n_stim_shuffle,'color',[0.5,0,0],'LineWidth',3)
title('Vm distribution residuals')

%%%%%%%%%%%%
SignalNorm = mean(DataNorm(1,:,:),2);
SignalShuffle = mean(DataShuffle(1,:,:),2);

edges = [0:0.5:20];
bytrial = squeeze(SignalNorm(1,:,1:sigon));
skew_pre_norm = skewness(bytrial(:))
n_pre_norm = histc(bytrial(:),edges) / (sigon*ntrials);
bytrial = squeeze(SignalNorm(1,:,sigon:sigoff));
skew_stim_norm = skewness(bytrial(:))
n_stim_norm = histc(bytrial(:),edges) / ((sigoff-sigon)*ntrials);

edges = [0:0.5:20];
bytrial = squeeze(SignalShuffle(1,:,1:sigon));
skew_pre_shuffle = skewness(bytrial(:))
n_pre_shuffle = histc(bytrial(:),edges) / (sigon*ntrials);
bytrial = squeeze(SignalShuffle(1,:,sigon:sigoff));
skew_stim_shuffle = skewness(bytrial(:))
n_stim_shuffle = histc(bytrial(:),edges) / ((sigoff-sigon)*ntrials);

figure;hold on
stairs(edges,n_pre_norm,'color',[0,0,0],'LineWidth',5)
stairs(edges,n_pre_shuffle,'color',[0.8,0.8,0.8],'LineWidth',3)
stairs(edges,n_stim_norm,'color',[1,0,0],'LineWidth',5)
stairs(edges,n_stim_shuffle,'color',[0.5,0,0],'LineWidth',3)
title('Vm distribution signal')

%%%%%%%%%%%
edges = [0:0.5:20];
bytrial = squeeze(DataNorm(1,:,1:sigon));
skew_pre_norm = skewness(bytrial(:))
n_pre_norm = histc(bytrial(:),edges) / (sigon*ntrials);
bytrial = squeeze(DataNorm(1,:,sigon:sigoff));
skew_stim_norm = skewness(bytrial(:))
n_stim_norm = histc(bytrial(:),edges) / ((sigoff-sigon)*ntrials);

edges = [0:0.5:20];
bytrial = squeeze(DataShuffle(1,:,1:sigon));
skew_pre_shuffle = skewness(bytrial(:))
n_pre_shuffle = histc(bytrial(:),edges) / (sigon*ntrials);
bytrial = squeeze(DataShuffle(1,:,sigon:sigoff));
skew_stim_shuffle = skewness(bytrial(:))
n_stim_shuffle = histc(bytrial(:),edges) / ((sigoff-sigon)*ntrials);

figure;hold on
stairs(edges,n_pre_norm,'color',[0,0,0],'LineWidth',5)
stairs(edges,n_pre_shuffle,'color',[0.8,0.8,0.8],'LineWidth',3)
stairs(edges,n_stim_norm,'color',[1,0,0],'LineWidth',5)
stairs(edges,n_stim_shuffle,'color',[0.5,0,0],'LineWidth',3)
title('Vm distribution total')

n_pre_shuffle = histc(TrialVarShuffle(1:sigon),edges) / sigon;
n_stim_shuffle = histc(TrialVarShuffle(sigon:sigoff),edges) / (sigoff-sigon);

TrialVarNorm = var(squeeze(DataNorm(1,:,:)),1);
TrialVarShuffle = var(squeeze(DataShuffle(1,:,:)),1);

figure;hold on
line(xtime_spikes,TrialVarNorm','color','k')
line(xtime_spikes,TrialVarShuffle','color','r')

sigon = min(find(xtime_spikes>0));
sigoff = max(find(xtime_spikes<6));

edges = [0:1:14];
n_pre_norm = histc(TrialVarNorm(1:sigon),edges) / sigon;
n_stim_norm = histc(TrialVarNorm(sigon:sigoff),edges) / (sigoff-sigon);

n_pre_shuffle = histc(TrialVarShuffle(1:sigon),edges) / sigon;
n_stim_shuffle = histc(TrialVarShuffle(sigon:sigoff),edges) / (sigoff-sigon);

figure;hold on
stairs(edges, n_pre_norm, 'color','k','LineWidth',5)
stairs(edges, n_stim_norm, 'color','r','LineWidth',3)
xlabel('trial var')
title('norm data; var prestim (black) and during stim (red)')

figure;hold on
stairs(edges, n_pre_shuffle, 'color','k','LineWidth',5)
stairs(edges, n_stim_shuffle, 'color','r','LineWidth',3)
xlabel('trial var')
title('shuffle data; var prestim (black) and during stim (red)')

figure;hold on
line(xtime_spikes,(TrialVarNorm./TrialVarShuffle)')
ylabel('TrialVarNorm./TrialVarShuffle')
xlabel('sec')

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

