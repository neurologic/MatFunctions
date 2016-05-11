cd('/Users/kperks/mnt/cube/Ice/kperks/B1087/')

% addpath('/mnt/cube/Ice/kperks/MatFunctions/');
% addpath('/home/kperks/Ice/MatFunctions/');

addpath('/Users/kperks/mnt/cube/Ice/kperks/MatFunctions/');


%% load up intracellular data
IntracellularData = load('IntracellularData.mat');
IntracellularData = IntracellularData.IntracellularData;

nsite_I = size(IntracellularData,2);

%% load up extracellular data
metatoes = load('metatoes.mat');

allsortclass = metatoes.allsortclass;
allchanxz = metatoes.allchanxz;
allclusterID = metatoes.allclusterid;
allsiteID = metatoes.allsiteID;


metatoes = metatoes.metatoes;

siteID_E = unique(allsiteID);
nsite_E = size(siteID_E,2);

siteID_E = unique(allsiteID);
nsite_E = size(siteID_E,2);

% get spikes mat filtered with synaptic
trial_edges = [-2,8.5];
exptdt = IntracellularData{1}.dt;
type = 'synaptic';
trialnsamps = round(diff(trial_edges)*(1/exptdt));
%  trialnsamps = 10000;
xtime_spikes = linspace(trial_edges(1),trial_edges(2),trialnsamps);
allSpikes = SpikeMat_filtered(metatoes,type);
% allSpikes is a cell array of function handles to pass xtime_spikes to for
% individual trials to build sub-matrices to analyze

fs = 31250;
stimfs = 40000;
clusterref = 1;
stimref = 1;
stimdur = (metatoes{clusterref}.stims{stimref}.stim_end_times(1)/fs)-(metatoes{clusterref}.stims{stimref}.stim_start_times(1)/fs);
stimstart = min(find(xtime_spikes >=0));
stimstop = min(find(xtime_spikes >=stimdur));


ncluster = size(allSpikes,1)
nstim = size(allSpikes,2)
ntrial = size(allSpikes,3)
trialnsamps

%% for a single stimulus, population response across trials across samples (for a simultaneously recorded population)
siteID = 3;
simulataneousInds = allSpikes(find(allsiteID == siteID),:,:);
stimID = 1;
for istim = stimID
    fs = 31250;
    stimfs = 40000;
    
    this_pop_response = zeros(ntrial,trialnsamps);
    for itrial = 1:ntrial
        this_trial_response = zeros(ncluster,trialnsamps);
        for icluster = 1:size(simulataneousInds,1)
            
            this_cluster_stim_trial_fun = simulataneousInds{icluster,istim,itrial};
            if ~isempty(this_cluster_stim_trial_fun)
                this_trial_response(icluster,:) = this_cluster_stim_trial_fun(xtime_spikes);
                
            end
            
        end
        this_pop_response(itrial,:) = sum(this_trial_response);
        
        this_pop_cluster_response(itrial,:,:) = (this_trial_response);
    end
    
end

%reshape by 50ms bins 
signal_pop_response = this_pop_response(:,(stimstart+diff([stimstart,stimstop])/6):stimstop);
numbins = 5000 / 100;
binsize = round(size(signal_pop_response,2)/numbins);
evenstart_ind = mod(size(signal_pop_response,2),1216) +1;
signal_pop_response = signal_pop_response(:,evenstart_ind:end);
numbins = size(signal_pop_response,2)/binsize;
pop_sum_response_by_bin = reshape(signal_pop_response,size(signal_pop_response,1),numbins,binsize);

trial1 = pop_sum_response_by_bin(1,:,:);

save('WorkingData/this_pop_response.mat','this_pop_response')
save('WorkingData/pop_sum_response_by_bin.mat','pop_sum_response_by_bin')

%reshape by 50ms bins 
signal_pop_cluster_response = this_pop_cluster_response(:,:,(stimstart+diff([stimstart,stimstop])/6):stimstop);
numbins = 5000 / 100;
binsize = round(size(signal_pop_cluster_response,2)/numbins);
evenstart_ind = mod(size(signal_pop_cluster_response,2),1216) +1;
signal_pop_cluster_response = signal_pop_cluster_response(:,evenstart_ind:end);
numbins = size(signal_pop_cluster_response,2)/binsize;
pop_sum_response_by_bin = reshape(signal_pop_cluster_response,size(signal_pop_cluster_response,1),numbins,binsize);

trial1 = pop_sum_response_by_bin(1,:,:);

save('WorkingData/this_pop_response.mat','this_pop_response')
save('WorkingData/pop_sum_response_by_bin.mat','pop_sum_response_by_bin')

%% for a single stimulus, intracellular response across trials across samples

cellID = 3;
thiscell = IntracellularData{cellID};

thiscell_thisstim = thiscell.stim_data{stimID


thisresponse = thiscell_thisstim(:,round(size(thiscell_thisstim,2)/6):end);
numbins = 5000 / 100;
binsize = round(size(thisresponse,2)/numbins);
evenstart_ind = mod(size(thisresponse,2),binsize) +1;
thisresponse = thisresponse(:,evenstart_ind:end);
numbins = size(thisresponse,2)/binsize;
intracell_response_by_bin = reshape(thisresponse,size(thisresponse,1),numbins,binsize);

trial1 = intracell_response_by_bin(1,:,:);

save('WorkingData/intracell_response_by_bin.mat','intracell_response_by_bin')
