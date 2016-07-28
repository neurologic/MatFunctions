cd('/mnt/cube/Ice/kperks/B1075/MetaMat/')
currentpath = pwd
% addpath('/mnt/cube/Ice/kperks/MatFunctions/');
%%%%%%%%%%%%%%%%
%load and set up data for analysis
%%%%%%%%%%%%%%%%

% path_to_phy = 'Z:\Ice\kperks\B1112\klusta\phy041316\';
path_to_phy = '../klusta/phy041216/';
sites = [];
d = dir(path_to_phy);
for id = 1:size(d,1)-2
    sites{id} = d(id+2).name;
end

%%%%%
% get indices for narrow and wide units
%%%%%
clusterID = [];
siteID = [];
isNarrow = [];
isNarrow_TP = [];
isNarrow_HH = [];
for id = 1:size(sites,2)
    load([path_to_phy sites{id} '/waveform_widths.mat'])
    thesesites = zeros(1,size(struct.cluster,2));
    thesesites(:) = id;
    siteID = [siteID, thesesites];
    clusterID = [clusterID, struct.cluster];
    isNarrow_TP = [isNarrow_TP, struct.narrow_trough_peak];
    isNarrow_HH = [isNarrow_HH, struct.narrow_half_height];

end
isNarrow = isNarrow_TP .* isNarrow_HH;

%%%%%
% load metatoes info
%%%%%
metatoes = load('metatoes.mat');
metatoesVars = metatoes;
vars = fieldnames(metatoes);
for ifield = 1:size(vars,1)
    s = [vars{ifield} '= metatoesVars.' vars{ifield} ';'];
    eval(s)
end

%%%%%
% calculate windows for response
%%%%%
fs = metatoes{1}.fs;
stimendtime = metatoes{1}.stims{1}.stim_end_times-metatoes{1}.stims{1}.stim_start_times;
stimdur = stimendtime(1)/fs;

sigon = min(find(xtime_spikes>0));
sigoff = max(find(xtime_spikes<stimdur));

calcstart = sigon + round((sigoff-sigon)/6);

%%%%%
% get indices for clusters want to include in analysis based on matatoes informatoin
%%%%%
% Good Quality
% From one penetration
qualinds = find(allsortclass == 2);
siteinds = find(allsiteID); %find(allsiteID == 4 | allsiteID == 5 | allsiteID == 6);

%%%%%%%%%%%%%%%%
% Wide Spike Clusters
%%%%%%%%%%%%%%%%
celltypeinds = find(isNarrow == 0);
include_inds = intersect(intersect(qualinds,siteinds),celltypeinds);

nclusters = size(include_inds,1);

nstims = size(allSpikes_SynFilter,2);
ntrials = size(allSpikes_SynFilter,3);
nsamps = size(xtime_spikes,2);

TheseSpikes = squeeze(allSpikes_SynFilter(include_inds,:,:));
these_locations = allchanxz(include_inds);

SignalFiltR = zeros(nclusters,nstims,sigoff-calcstart+1);
%get trial-averaged response per clsuter per stimulus
for isite = 1:nclusters
    siteprogress_countdown = nclusters - isite;

    fileID = fopen('matlog.txt','w');
    t = datestr(datetime('now'));
    fprintf(fileID,'%s %i %s\n','Wide_countdown_Prep',siteprogress_countdown,t);
    fclose(fileID);
    
        for istim = 1:nstims
          
            TrialFiltR = zeros(ntrials,nsamps);
            for itrial = 1:ntrials
                trialfun1 = TheseSpikes{isite,istim,itrial};
                if ~isempty(trialfun1)
                    TrialFiltR(itrial,:) = trialfun1(xtime_spikes); %for this cluster on this trial
                end
            end
            SignalFiltR(isite,istim,:) = mean(TrialFiltR(:,calcstart:sigoff));
        end
    
end
save('TrialAveragedResponse_WideSpikes.mat','SignalFiltR','include_inds')

SignalCorr = nan(nclusters,nclusters,nstims);
PairWiseDistance = nan(nclusters,nclusters);

%need specific index for corr matrix because sometimes all NaN
idx = true(2);
idx = ~(tril(idx));

for isite1 = 1:nclusters
    siteprogress_countdown = nclusters - isite1;

    fileID = fopen('matlog.txt','w');
    t = datestr(datetime('now'));
    fprintf(fileID,'%s %i %s\n','Wide_countdown',siteprogress_countdown,t);
    fclose(fileID);

    for isite2 = isite1:nclusters
        
        for istim = 1:nstims
          
%             TrialFiltR = zeros(2,ntrials,nsamps);
%             for itrial = 1:ntrials
%                 trialfun1 = TheseSpikes{isite1,istim,itrial};
%                 if ~isempty(trialfun1)
%                     TrialFiltR(1,itrial,:) = trialfun1(xtime_spikes); %for this cluster on this trial
%                 end
%                  trialfun2 = TheseSpikes{isite2,istim,itrial};
%                 if ~isempty(trialfun2)
%                     TrialFiltR(2,itrial,:) = trialfun2(xtime_spikes); %for this cluster on this trial
%                 end
%                 
%             end
            c = corrcoef(squeeze(SignalFiltR(isite1,istim,:)),squeeze(SignalFiltR(isite2,istim,:)));
           
            SignalCorr(isite1,isite2,istim) = c(idx);
            
            
        end
        PairWiseDistance(isite1,isite2) = abs(these_locations(isite2) - these_locations(isite1));
    end
end

save('SpatialSignalCorr_WideSpikes.mat','PairWiseDistance','SignalCorr','include_inds')

%%%%%%%%%%%%%%%%
% Narrow Spike Clusters
%%%%%%%%%%%%%%%%
celltypeinds = find(isNarrow == 1);
include_inds = intersect(intersect(qualinds,siteinds),celltypeinds);

nclusters = size(include_inds,1);

nstims = size(allSpikes_SynFilter,2);
ntrials = size(allSpikes_SynFilter,3);
nsamps = size(xtime_spikes,2);

TheseSpikes = squeeze(allSpikes_SynFilter(include_inds,:,:));
these_locations = allchanxz(include_inds);

SignalFiltR = zeros(nclusters,nstims,sigoff-calcstart+1);
%get trial-averaged response per clsuter per stimulus
for isite = 1:nclusters
    
    siteprogress_countdown = nclusters - isite;

    fileID = fopen('matlog.txt','w');
    t = datestr(datetime('now'));
    fprintf(fileID,'%s %i %s\n','Narrow_countdown prep',siteprogress_countdown,t);
    fclose(fileID);
    
        for istim = 1:nstims
          
            TrialFiltR = zeros(ntrials,nsamps);
            for itrial = 1:ntrials
                trialfun1 = TheseSpikes{isite,istim,itrial};
                if ~isempty(trialfun1)
                    TrialFiltR(itrial,:) = trialfun1(xtime_spikes); %for this cluster on this trial
                end
            end
            SignalFiltR(isite,istim,:) = mean(squeeze(TrialFiltR(:,calcstart:sigoff)));
        end
    
end
save('TrialAveragedResponse_NarrowSpikes.mat','SignalFiltR','include_inds')

SignalCorr = nan(nclusters,nclusters,nstims);
PairWiseDistance = nan(nclusters,nclusters);

%need specific index for corr matrix because sometimes all NaN
idx = true(2);
idx = ~(tril(idx));

for isite1 = 1:nclusters
    siteprogress_countdown = nclusters - isite1;

    fileID = fopen('matlog.txt','w');
    t = datestr(datetime('now'));
    fprintf(fileID,'%s %i %s\n','Narrow_countdown',siteprogress_countdown,t);
    fclose(fileID);

    for isite2 = isite1:nclusters
        
        for istim = 1:nstims
          
%             TrialFiltR = zeros(2,ntrials,nsamps);
%             for itrial = 1:ntrials
%                 trialfun1 = TheseSpikes{isite1,istim,itrial};
%                 if ~isempty(trialfun1)
%                     TrialFiltR(1,itrial,:) = trialfun1(xtime_spikes); %for this cluster on this trial
%                 end
%                  trialfun2 = TheseSpikes{isite2,istim,itrial};
%                 if ~isempty(trialfun2)
%                     TrialFiltR(2,itrial,:) = trialfun2(xtime_spikes); %for this cluster on this trial
%                 end
%                 
%             end
            c = corrcoef(squeeze(SignalFiltR(isite1,istim,:)),squeeze(SignalFiltR(isite2,istim,:)));
           
            SignalCorr(isite1,isite2,istim) = c(idx);
            
            
        end
        PairWiseDistance(isite1,isite2) = abs(these_locations(isite2) - these_locations(isite1));
    end
end

save('SpatialSignalCorr_NarrowSpikes.mat','PairWiseDistance','SignalCorr','include_inds')


