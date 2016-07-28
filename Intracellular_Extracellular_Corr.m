% with an intracellular cell and an extracellular dataset

% 1 - bin trial-averaged response (take mean Vm per bin)
% across all locations for Vm (intracellular response)
% at each location for PVm (pseudo intracellular - clusters within Z microns of location)
% 2 -

% correlate spiking activity of intracellular cell with spiking activity of
% all clusters and find the locations with correlation in top 1% or
% something (... to correlate with predicted corr for close
% is this consistent with where Vm is mose closely correlated with?
cd('/mnt/cube/Ice/kperks/B1087/MetaMat/')
currentpath = pwd

path_to_phy = '../klusta/phy040516/';

CellID = 1; %which intracellular cell

%%%%%
% Get extracellular pop sliding window
%%%%%
DeltaD = 100;
BreadthD = 100;

%%%%%
% Load_Metatoes
%%%%%
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
IntracellStimNames = IntracellularData{CellID}.stimnames;

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
siteinds = find(allsiteID); %find(allsiteID == 1 | allsiteID == 2 | allsiteID == 3);

%%%%%
% Set up locations and window bins
%%%%%
DCenters = [min(unique(allchanxz))+BreadthD:DeltaD:max(unique(allchanxz))-BreadthD];
ncenters = size(DCenters,2);

LocalPop = zeros(ncenters,nstims,sigoff-sigon+1);
% LocalPop_Wide = zeros(ncenters,nstims,sigoff-sigon+1);
% LocalPop_Narrow = zeros(ncenters,nstims,sigoff-sigon+1);
for icenter = 1:ncenters
    fileID = fopen('matlog.txt','w');
    t = datestr(datetime('now'));
    fprintf(fileID,'%s %i %s\n','calculating location',icenter,t);
    fclose(fileID);
    
    locationinds = intersect(find(allchanxz>=(DCenters(icenter)-BreadthD)),find(allchanxz<=(DCenters(icenter)+BreadthD)));
    
    %%%%%%%%%%%%%%%%
    % Wide Spike Clusters
    %%%%%%%%%%%%%%%%
    celltypeinds = find(isNarrow == 0);
    include_inds_wide = intersect(intersect(intersect(qualinds,siteinds),celltypeinds),locationinds);
    
    nclusters = size(include_inds_wide,1);
    
    nstims = size(allSpikes_SynFilter,2);
    ntrials = size(allSpikes_SynFilter,3);
    nsamps = size(xtime_spikes,2);
    
    TheseSpikes = (allSpikes_SynFilter(include_inds_wide,:,:));
    these_locations = allchanxz(include_inds_wide);
    
    SignalFiltR_wide = zeros(nclusters,nstims,sigoff-sigon+1);
    %get trial-averaged response per clsuter per stimulus
    for isite = 1:nclusters
        
        for istim = 1:nstims
            
            TrialFiltR = zeros(ntrials,nsamps);
            for itrial = 1:ntrials
                trialfun1 = TheseSpikes{isite,istim,itrial};
                if ~isempty(trialfun1)
                    TrialFiltR(itrial,:) = trialfun1(xtime_spikes); %for this cluster on this trial
                end
            end
            SignalFiltR_wide(isite,istim,:) = mean(TrialFiltR(:,sigon:sigoff));
        end
        
    end
    
%     LocalPop_Wide(icenter,:,:) = squeeze(mean(SignalFiltR_wide,1));
    
    %%%%%%%%%%%%%%%%
    % Narrow Spike Clusters
    %%%%%%%%%%%%%%%%
    celltypeinds = find(isNarrow == 1);
    include_inds_narrow = intersect(intersect(intersect(qualinds,siteinds),celltypeinds),locationinds);
    
    nclusters = size(include_inds_narrow,1);
    
    TheseSpikes = (allSpikes_SynFilter(include_inds_narrow,:,:));
    these_locations = allchanxz(include_inds_narrow);
    
    SignalFiltR_narrow = zeros(nclusters,nstims,sigoff-sigon+1);
    %get trial-averaged response per clsuter per stimulus
    for isite = 1:nclusters
        for istim = 1:nstims
            
            TrialFiltR = zeros(ntrials,nsamps);
            for itrial = 1:ntrials
                trialfun1 = TheseSpikes{isite,istim,itrial};
                if ~isempty(trialfun1)
                    TrialFiltR(itrial,:) = trialfun1(xtime_spikes); %for this cluster on this trial
                end
            end
            SignalFiltR_narrow(isite,istim,:) = mean(TrialFiltR(:,sigon:sigoff),1);
        end
    end
    
    popsum = sum([(SignalFiltR_wide(:,:,:));-(SignalFiltR_narrow(:,:,:))],1);
    LocalPop(icenter,:,:) = squeeze(popsum);
end

save(['IntraExtraCorrLags_Delta' num2str(DeltaD) '_Breadth' num2str(BreadthD) '.mat'],'LocalPop','DCenters','BreadthD','ExtracellStimNames','IntracellStimNames')
%%
CellID = 2;
ncenters = size(DCenters,2);
nstims = size(LocalPop,2);
medsmoothing = 400;
% LocationLag_narrow = nan(ncenters,nstims);
LocationLag = nan(ncenters,nstims);
% LocationCorr_narrow = nan(ncenters,nstims);
LocationCorr = nan(ncenters,nstims);
% LocationCorr_WideNarrow = nan(ncenters,nstims);
% LocationCorr_WideNarrow  = nan(ncenters,nstims);
for icenter = 1:ncenters
    
    for istim = 1:nstims
        Vm = mean(IntracellularData{CellID}.stim_data{istim},1);
        Vm = Vm - min(Vm);
        Vm = medfilt1(Vm,medsmoothing,[],2);
        Vm = Vm(1:size(LocalPop,3));
        
        PVm = medfilt1(squeeze(LocalPop(icenter,istim,:)),medsmoothing,[],1);
%         PVm_narrow = medfilt1(squeeze(LocalPop_Narrow(icenter,istim,:)),medsmoothing,[],1);
        
        [r,lags] = xcorr(Vm,PVm,'coeff');
        maxind = find(r == max(r));
        if ~isempty(maxind)
            LocationCorr(icenter,istim) = r(maxind);
            LocationLag(icenter,istim) = lags(maxind);
        end
        
%         [r,lags] = xcorr(Vm,PVm_narrow,'coeff');
%         maxind = find(r == max(r));
%         if ~isempty(maxind)
%             LocationCorr_narrow(icenter,istim) = r(maxind);
%             LocationLag_narrow(icenter,istim) = lags(maxind);
%         end
        
%         [r,lags] = xcorr(PVm_wide,PVm_narrow,'coeff');
%         maxind = find(r == max(r));
%         if ~isempty(maxind)
%             LocationCorr_WideNarrow(icenter,istim) = r(maxind);
%             LocationLag_WideNarrow(icenter,istim) = lags(maxind);
%         end
    end
    
end


ncenters = size(DCenters,2);
figure;
hold on
for icenter = 1:ncenters
scatter(repmat(DCenters(icenter),1,size(LocationCorr,2)),LocationCorr(icenter,:),100,'k','fill')
scatter(DCenters(icenter),mean(LocationCorr(icenter,:)),400,'r','fill')
% scatter(repmat(DCenters(icenter),1,size(LocationCorr_narrow,2)),LocationCorr_narrow(icenter,:),100,'r','fill')
end
title(IntracellularData{CellID}.exptname,'Interpreter','none')
xlabel('D binned at each 100 with 200u width')
ylabel('max corr')


ncenters = size(DCenters,2);
figure;
hold on
for icenter = 1:ncenters
scatter(repmat(DCenters(icenter),1,size(LocationLag,2)),LocationLag(icenter,:),100,'k','fill')
scatter(DCenters(icenter),mean(LocationLag(icenter,:)),400,'r','fill')
end
title(IntracellularData{CellID}.exptname,'Interpreter','none')
xlabel('D binned at each 100 with 200u width')
ylabel('lag at max corr')

figure;scatter(LocationCorr(:),LocationLag(:),100,'k','fill')
xlabel('All Max Corr')
ylabel('All Lags at max Corr')



% save(['IntraExtraCorrLags_Delta' num2str(DeltaD) '_Breadth' num2str(BreadthD) '.mat'],'LocationLag','LocationCorr','LocalPop','DCenters','BreadthD','ExtracellStimNames','IntracellStimNames')
%%
% 
IECorrs = load('IntraExtraCorrLags_Delta100_Breadth100.mat');
vars = fieldnames(IECorrs);
for ifield = 1:size(vars,1)
    s = [vars{ifield} '= IECorrs.' vars{ifield} ';'];
    eval(s)
end

ncenters = size(DCenters,2);
figure;
hold on
for icenter = 1:ncenters
scatter(repmat(DCenters(icenter),1,size(LocationCorr,2)),LocationCorr(icenter,:),100,'k','fill')
% scatter(repmat(DCenters(icenter),1,size(LocationCorr_narrow,2)),LocationCorr_narrow(icenter,:),100,'r','fill')
end

figure;
hold on
for icenter = 1:ncenters
scatter(repmat(DCenters(icenter),1,size(LocationLag_wide,2)),LocationLag_wide(icenter,:),100,'b','fill')
scatter(repmat(DCenters(icenter),1,size(LocationLag_narrow,2)),LocationLag_narrow(icenter,:),100,'r','fill')
end

figure;
hold on
for icenter = 1:ncenters
scatter(repmat(DCenters(icenter),1,size(LocationCorr_WideNarrow,2)),LocationCorr_WideNarrow(icenter,:),100,'k','fill')
end
figure;
hold on
for icenter = 1:ncenters
scatter(repmat(DCenters(icenter),1,size(LocationLag_WideNarrow,2)),LocationLag_WideNarrow(icenter,:),100,'k','fill')
end
%%

figure; 
plot(sum([squeeze(SignalFiltR_wide(:,1,:));-squeeze(SignalFiltR_narrow(:,1,:))],1)')
load('IntraExtraCorrLags.mat','LocationLag_wide','LocationCorr_wide','LocationLag_narrow','LocationCorr_narrow','LocalPop_Narrow','LocalPop_Wide','DCenters','BreadthD','ExtracellStimNames')