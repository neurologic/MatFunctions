cd('/mnt/cube/Ice/kperks/B1087/MetaMat/')
currentpath = pwd

path_to_phy = '../klusta/phy040516/';
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

include_inds = intersect(qualinds,siteinds);
TheseNarrowIDs = isNarrow(include_inds);
TheseSites = allsiteID(include_inds);

nclusters = size(include_inds,1);

nstims = size(allSpikes_SynFilter,2);
ntrials = size(allSpikes_SynFilter,3);
nsamps = size(xtime_spikes,2);

TheseSpikes = squeeze(allSpikes_SynFilter(include_inds,:,:));
these_locations = allchanxz(include_inds);

SignalFiltR = zeros(nclusters,nstims,nsamps);
%get trial-averaged response per clsuter per stimulus
for istim = 1;%2:nstims
    TrialFiltR = zeros(nclusters,ntrials,nsamps);
    
    for icluster = 1:nclusters
        
        fileID = fopen('matlog.txt','w');
        t = datestr(datetime('now'));
        fprintf(fileID,'%s %i %s %i %s\n','istim',istim,'isite',icluster,t);
        fclose(fileID);
        
        doflip = 0;
        if TheseNarrowIDs(icluster)==1
            doflip = 1;
        end
        
        for itrial = 1:ntrials
            trialfun1 = TheseSpikes{icluster,istim,itrial};
            if ~isempty(trialfun1)
                TrialFiltR(icluster,itrial,:) = trialfun1(xtime_spikes); %for this cluster on this trial
            end
        end
    end
    
    eachsite = unique(TheseSites);
    
    SiteResp = nan(size(eachsite,2),ntrials,nsamps);
    for isite = 1:size(eachsite,2)
        theseinds = find(TheseSites==eachsite(isite));
        SiteResp(isite,:,:) = (sum(TrialFiltR(theseinds,:,:),1));
        
    end
    
    s = ['TrialResponses_stim_' num2str(istim) '_wav_' metatoes{1}.stims{istim}.name{1}(1) '.mat'];
    save(['TrialResponsesByStim/' s],'SiteResp')
end



%%

% s = ['TrialResponses_stim_' num2str(istim) '_wav_' metatoes{1}.stims{istim}.name{1}(1) '.mat'];
%     save(['TrialResponsesByStim/' s],'TrialFiltR','include_inds','TheseNarrowIDs','TheseSites','-v7.3')

%%% write code to save each of the stimuli as a mat variable to be called
%%% individually that has the individual trial responses for that stimulus
%%% for each cluster

% save('TrialAveraged_WholeTrial.mat','SignalFiltR','NarrowInds','include_inds','sigon','sigoff','calcstart')