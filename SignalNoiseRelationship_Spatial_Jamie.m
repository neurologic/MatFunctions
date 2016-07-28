% for each site in each bird, do Jamie's Signal Noise Correlation Analysis
% start with 5 different data set segmentations:
    % 1 - by sample
    % 2 - by song (excluding first motif)
    % 3 - by 500msec
    % 4 - by 1sec
    % 5 - only the 500msecond after common intro
% if see difference with training then the first motif can be a control?
% analyze how this relationship is spatially dependent 

%%%%%
% parameters to determine before running analysis
%%%%%
response_bin_dur = 0.5; %seconds
siteID = 2;

%%%%% 
% Load_Metatoes
%%%%%

metatoes = load('metatoes.mat');
metatoesVars = metatoes;
vars = fieldnames(metatoes);
for ifield = 1:size(vars,1)
    s = [vars{ifield} '= metatoesVars.' vars{ifield} ';'];
    eval(s)
end

%%%%% 
% Load Behavior Stim data
%%%%%
load('behavior.mat','metatoes_songlist','trainedSong_ind','novelSong_ind')

%%%%%
% Decide which clusters from metatoes to include
%%%%%

these_classinds = find(allsortclass == 2);
these_siteinds = find(allsiteID == siteID);
include_inds = intersect(these_classinds,these_siteinds);
thesetoes = metatoes(include_inds);
theselocations = allchanxz(include_inds);

%%%%%
% calc params from data like ntrials and nclusters, etc
%%%%%
nclusters = size(include_inds,1);
nstims = size(metatoes{1}.stims,1);
ntrials = metatoes{1}.stims{1}.ntrials;

%%%%%
% get response windows and bin the response for each cluster
%%%%%
fs = metatoes{1}.fs;
stimendtime = metatoes{1}.stims{1}.stim_end_times-metatoes{1}.stims{1}.stim_start_times;
stimdur = stimendtime(1)/fs;
calcstart = round(stimdur/6);
% sigon = min(find(xtime_spikes>0));
% sigoff = max(find(xtime_spikes<stimdur));
% 
% calcstart = sigon + round((sigoff-sigon)/6);

bins = [calcstart:response_bin_dur:stimdur];
nbins = size(bins,2);
%%%%%
% build binned response mat (prep) 
%%%%%

BinnedSpikeR = zeros(nclusters,nstims,nbins-1,ntrials);
for icluster = 1:nclusters

    for isong = 1:nstims

        for itrial = 1:ntrials
            thistrial = thesetoes{icluster}.stims{isong}.toes{itrial};
            
            for ibin = 1:nbins-1
                numspikes = 0;
                spikeinds = intersect(find(thistrial>=bins(ibin)),find(thistrial<bins(ibin+1)));
                if ~isempty(spikeinds)
                    numspikes = size(spikeinds,1);
                end
                BinnedSpikeR(icluster,isong,ibin,itrial) = numspikes/response_bin_dur;
            end
            
        end
    end
end


idx = true(2);
idx = ~tril(idx);

%trained songs
nstims = size(trainedSong_ind,2);

TrialCorr_Trained = nan(nclusters,nclusters,(nbins-1)*nstims);
SignalHat_Trained = zeros(nclusters,(nbins-1)*nstims);
PairWiseDistance = nan(nclusters,nclusters);
for isite1 = 1:nclusters
    siteprogress_countdown = nclusters - isite1

    for isite2 = isite1:nclusters
        
        stimind = 1;
        for isong = 1:nstims
            for ibin = 1:nbins-1
                R1 = squeeze(BinnedSpikeR(isite1,trainedSong_ind(isong),ibin,:));
                R2 = squeeze(BinnedSpikeR(isite2,trainedSong_ind(isong),ibin,:));
                
                cT = corrcoef(R1,R2);
                TrialCorr_Trained(isite1,isite2,stimind) = cT(idx);
                
                SignalHat_Trained(isite1,stimind) = nanmean(R1);
                
                stimind = stimind +1;
            end
        end
        PairWiseDistance(isite1,isite2) = abs(theselocations(isite2) - theselocations(isite1));
    end
end
SignalCorr_Trained = nan(nclusters,nclusters);
for isite1 = 1:nclusters
    siteprogress_countdown = nclusters - isite1

    for isite2 = isite1:nclusters
        cS = corrcoef(SignalHat_Trained(isite1,:),SignalHat_Trained(isite2,:));
        SignalCorr_Trained(isite1,isite2) = cS(idx);
        
    end
end
TrialCorrHat_Trained = squeeze(nanmean(TrialCorr_Trained,3));

idx = true(nclusters);
idx = ~tril(idx);

figure;hold on
scatter(SignalCorr_Trained(idx),TrialCorrHat_Trained(idx))
title('trained stims')
xlabel('mean Signal Corr per pair')
ylabel('mean Trial Corr per pair')
SC_Trained = SignalCorr_Trained(idx);
TC_Trained = TrialCorrHat_Trained(idx);
cST_trained = corrcoef(SC_Trained(~isnan(TC_Trained)),TC_Trained(~isnan(TC_Trained)))
set(gca,'XLim',[-0.5,1])
set(gca,'YLim',[-0.5,1])

idx = true(2);
idx = ~tril(idx);
%novel songs
nstims = size(novelSong_ind,2);
TrialCorr_Novel = nan(nclusters,nclusters,(nbins-1)*nstims);
SignalHat_Novel = zeros(nclusters,(nbins-1)*nstims);
% PairWiseDistance = nan(nclusters,nclusters);
for isite1 = 1:nclusters
    siteprogress_countdown = nclusters - isite1

    for isite2 = isite1:nclusters
        
        stimind = 1;
        for isong = 1:nstims
            for ibin = 1:nbins-1
                R1 = squeeze(BinnedSpikeR(isite1,novelSong_ind(isong),ibin,:));
                R2 = squeeze(BinnedSpikeR(isite2,novelSong_ind(isong),ibin,:));
                
                cT = corrcoef(R1,R2);
                TrialCorr_Novel(isite1,isite2,stimind) = cT(idx);
                
                SignalHat_Novel(isite1,stimind) = nanmean(R1);
                
                stimind = stimind +1;
            end
        end
%         PairWiseDistance(isite1,isite2) = abs(theselocations(isite2) - theselocations(isite1));
    end
end
SignalCorr_Novel = nan(nclusters,nclusters);
for isite1 = 1:nclusters
    siteprogress_countdown = nclusters - isite1

    for isite2 = isite1:nclusters
        cS = corrcoef(SignalHat_Novel(isite1,:),SignalHat_Novel(isite2,:));
        SignalCorr_Novel(isite1,isite2) = cS(idx);
        
    end
end
TrialCorrHat_Novel = squeeze(nanmean(TrialCorr_Novel,3));

idx = true(nclusters);
idx = ~tril(idx);

figure;hold on
scatter(SignalCorr_Novel(idx),TrialCorrHat_Novel(idx))
title('novel stims')
xlabel('mean Signal Corr per pair')
ylabel('mean Trial Corr per pair')
SC_Novel = SignalCorr_Novel(idx);
TC_Novel = TrialCorrHat_Novel(idx);
cST_novel = corrcoef(SC_Novel(~isnan(TC_Novel)),TC_Novel(~isnan(TC_Novel)))
set(gca,'XLim',[-0.5,1])
set(gca,'YLim',[-0.5,1])

corrEdges = [-1:0.01:1];
n_S_T = histc(SC_Trained,corrEdges)./size(SC_Trained,1);
n_S_N = histc(SC_Novel,corrEdges)./size(SC_Novel,1);
figure;hold on
stairs(corrEdges,n_S_T,'color','r')
stairs(corrEdges,n_S_N,'color','b')
title('signal correlation')
legend('trained','novel')

n_T_T = histc(TC_Trained,corrEdges)./size(TC_Trained,1);
n_T_N = histc(TC_Novel,corrEdges)./size(TC_Novel,1);
figure;hold on
stairs(corrEdges,n_T_T,'color','r')
stairs(corrEdges,n_T_N,'color','b')
title('trial correlation')
legend('trained','novel')

REdges = [0:0.5:ceil(max([max(max(SignalHat_Novel)),max(max(SignalHat_Trained))]))];
n_T = histc(SignalHat_Trained(:),REdges)./size(SignalHat_Trained(:),1);
n_N = histc(SignalHat_Novel(:),REdges)./size(SignalHat_Novel(:),1);
figure;hold on
stairs(REdges,n_T,'color','r')
stairs(REdges,n_N,'color','b')
title('Mean Signal Response (each bin each cluster)')
legend('trained','novel')
axis tight

REdges = [1:0.5:ceil(max([max(max(SignalHat_Novel)),max(max(SignalHat_Trained))]))];
n_T = histc(SignalHat_Trained(:),REdges)./size(SignalHat_Trained(:),1);
n_N = histc(SignalHat_Novel(:),REdges)./size(SignalHat_Novel(:),1);
figure;hold on
stairs(REdges,n_T,'color','r')
stairs(REdges,n_N,'color','b')
title('Mean Signal Response (each bin each cluster)')
legend('trained','novel')
axis tight