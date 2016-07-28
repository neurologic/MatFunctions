cd('/mnt/cube/Ice/kperks/B1083/MetaMat/')
currentpath = pwd

path_to_phy = '../klusta/phy040516/';

% CellID = 1; %which intracellular cell
%%
load('IntracellularData.mat','IntracellularData')

metatoes = load('metatoes.mat');
metatoesVars = metatoes;
vars = fieldnames(metatoes);
for ifield = 1:size(vars,1)
    s = [vars{ifield} '= metatoesVars.' vars{ifield} ';'];
    eval(s)
end

for istim = 1:size(metatoes{1}.stims,1)
    ExtracellStimNames(istim) = metatoes{1}.stims{istim}.name;
end
IntracellStimNames = IntracellularData{1}.stimnames;

fs = metatoes{1}.fs;
stimendtime = metatoes{1}.stims{1}.stim_end_times-metatoes{1}.stims{1}.stim_start_times;
stimdur = stimendtime(1)/fs;
sigon = min(find(xtime_spikes>0));
sigoff = max(find(xtime_spikes<stimdur));

%%%%%
% get wide/narrow information
%%%%%
sites = [];
d = dir(path_to_phy);
for id = 1:size(d,1)-2
    sites{id} = d(id+2).name;
end

clusterID = [];
siteID = [];
isNarrow = [];
isNarrow_TP = [];
isNarrow_HH = [];
for id = 1:size(sites,2)
    waveparams = load([path_to_phy sites{id} '/waveform_widths.mat']);
    waveparams = waveparams.struct;
    thesesites = zeros(1,size(waveparams.cluster,2));
    thesesites(:) = id;
    siteID = [siteID, thesesites];
    clusterID = [clusterID, waveparams.cluster];
    isNarrow_TP = [isNarrow_TP, waveparams.narrow_trough_peak];
    isNarrow_HH = [isNarrow_HH, waveparams.narrow_half_height];
    
end
isNarrow = isNarrow_TP .* isNarrow_HH;

% Good Quality
% From all sites
qualinds = find(allsortclass == 2);
siteinds = find(allsiteID); %find(allsiteID == 1 | allsiteID == 2 | allsiteID == 3); %find(allsiteID); %

% celltypeinds = find(isNarrow == 0);
include_inds = intersect(qualinds,siteinds);

nclusters = size(include_inds,1);

nstims = size(allSpikes_SynFilter,2);
ntrials = size(allSpikes_SynFilter,3);
nsamps = size(xtime_spikes,2);

TheseSpikes = (allSpikes_SynFilter(include_inds,:,:));

%% if trial shuffling



type = 'synaptic';
% allSpikes_SynFilter = SpikeMat_filtered(metatoes,type);

toedata = metatoes(include_inds);
allSpikes = cell(size(toedata,1),size(toedata{1}.stims,1),size(toedata{1}.stims{1}.toes,1));
tau = 0.01;

alltoes = [];
ntoes = nan(nclusters,nstims,ntrials);
for isite = 1:size(toedata,1)
    
    thissite = toedata{isite};
    
    for istim = 1:max(size(toedata{isite}.stims));
        thisstim = thissite.stims{istim};
        
        for itrial = 1:thisstim.ntrials
            
            thisstim_toes = thisstim.toes;
            
            thistrial_toes = thisstim_toes{itrial};
            
            alltoes = [alltoes;thistrial_toes];
            
            ntoes(isite,istim,itrial) = size(thistrial_toes,1);
%             filteredspikes = zeros(1,nsamps);
        end
    end
end

shufflevec = randperm(size(alltoes,1));
shuffletoes = alltoes(shufflevec);

for isite = 1:size(toedata,1)
    
    
    for istim = 1:max(size(toedata{isite}.stims));
        
        
        for itrial = 1:thisstim.ntrials
            
            N = ntoes(isite,istim,itrial);
            
            thistrial_toes = [];
            if N > 0
                thistrial_toes = shuffletoes(1:N);
                shuffletoes = shuffletoes(N+1:end);
                
                if ~isempty(thistrial_toes)
                    
                    fun = filtered_response_synaptic(thistrial_toes, tau);
                    allSpikes{isite,istim,itrial,:} = fun;
                    
                    %                 filteredspikes = fun(linspace(trial_edges(1),trial_edges(2),nsamps));
                end
                
            end
            %             filteredspikes = zeros(1,nsamps);
        end
    end
end


TheseSpikes = allSpikes;

SignalFiltR = zeros(nclusters,nstims,sigoff-sigon+1);

%get trial-averaged response per clsuter per stimulus
for isite = 1:20% nclusters
    fileID = fopen('matlog.txt','w');
    t = datestr(datetime('now'));
    fprintf(fileID,'%s %i %s\n','calculating cluster',isite,t);
    fclose(fileID);
    for istim = 1:nstims
        
        TrialFiltR = zeros(ntrials,nsamps);
        for itrial = 1:ntrials
            trialfun1 = TheseSpikes{isite,istim,itrial};
            if ~isempty(trialfun1)
                TrialFiltR(itrial,:) = trialfun1(xtime_spikes); %for this cluster on this trial
            end
        end
        SignalFiltR(isite,istim,:) = mean(TrialFiltR(:,sigon:sigoff));
    end
    
end

save('MVRegressPrep_Shuffled.mat','SignalFiltR','include_inds','isNarrow','ExtracellStimNames','IntracellStimNames')
%%
% 
% Pop =[(SignalFiltR);-(SignalFiltR_narrow(:,1,:))];
% 
% 
% for istim = 1:nstims
%     Vm = mean(IntracellularData{CellID}.stim_data{istim},1);
%     Vm = Vm - min(Vm);
% %     Vm = medfilt1(Vm,medsmoothing,[],2);
%     Vm = Vm(1:size(Pop,3));
%     Vm = Vm';
%     
%     PVm = squeeze(Pop(:,istim,:));
%     PVm = PVm';
%     
%     [B,Sigma] = mvregress(PVm,Vm);
%     Beta(istim,:) = B;
% end
% 
% VmHat = [];
% for iterm = 1:size(B,1)
% VmHat(:,iterm) = B(iterm)*PVm(:,iterm);
% end

