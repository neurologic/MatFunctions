%% Given any SignalCorr and Distance matrices, get relationship
%matrices are [ncluster x ncluster x nstims]

%%%%TODO
% for each cluster pair for each stimulus, get the joint mean activity
% across the stimulus being considered and plot the ratio of signal
% corrlation to response magnitude. 

close all
do_visualize=1;
do_fit = 1;
do_CorrHat = 1;


BirdID = 'B1083';

load('behavior.mat')
songinds = union(trainedSong_ind,novelSong_ind);

matlist = {'SpatialSignalCorr_WideSpikes.mat','SpatialSignalCorr_NarrowSpikes.mat'};
for imat =1:2
    
    load(matlist{imat}, 'PairWiseDistance','SignalCorr','include_inds');
    
    
    % not all SignalCorr are real since some clusters do not respond to all
    MeanRs = squeeze(nanmean(SignalCorr(:,:,songinds),3));
    %get indices of upper triangle of corr matrix
    idx = true(size(MeanRs,1));
    idx = ~(tril(idx));
    % get vector for upper triangle of SignalCorr and Distance
    R = MeanRs(idx);
    D = PairWiseDistance(idx);
    
    %for each distance
    distancebins = unique(D);
    % will want to plot in a different color
    cols = colormap(hsv(size(distancebins,1)));
    
    %create figure to plot
    Rmean_bydistance = [];
    
    for iD = 1:size(distancebins,1)
        %each time this distance between clusters appears
        theseinds = find(D == distancebins(iD));
        Rmean_bydistance(iD) = nanmean(R(theseinds));
    end
    
    if imat ==1 % (so for wide)
        R_wide = R;
        D_wide = D;
        distancebins_wide = distancebins;
        Rmean_bydistance_wide = Rmean_bydistance;
    end
    if imat ==2 % (so for narrow)
        R_narrow = R;
        D_narrow = D;
        distancebins_narrow = distancebins;
        Rmean_bydistance_narrow = Rmean_bydistance;
    end
end

%%%%%%%%%%%%%%%%%
%eventually will not need to plot, but will plot now and look at what the
%metric of spatial dependence of signal correlation might be
%keep this code around to look manually at relationships in the future
%%%%%%%%%%%%%%%%%
%

if do_CorrHat ==1
    distance_width = 200;
    bininds_wide = (find(distancebins_wide<=distance_width));
    Rhat_wide = mean(Rmean_bydistance_wide(bininds_wide))
    
    bininds_narrow = (find(distancebins_narrow<=distance_width));
    Rhat_narrow = mean(Rmean_bydistance_narrow(bininds_narrow))
end
    
if do_visualize ==1
    figure;hold on
    scatter(D_wide,R_wide,10,[0.5 0.5 0.5])
    scatter(distancebins_wide,Rmean_bydistance_wide,100,'k','fill')
    axis tight
    set(gca,'YLim',[-0.2,1])
    title([BirdID '; Wide Spikes'])
    xlabel('distance between cluster')
    ylabel('mean signal correlation across songs')
    
    figure;hold on
    scatter(D_narrow,R_narrow,10,[0.5 0.5 0.5])
    scatter(distancebins_narrow,Rmean_bydistance_narrow,100,'k','fill')
    axis tight
    set(gca,'YLim',[-0.2,1])
    title([BirdID '; Narrow Spikes'])
    xlabel('distance between cluster')
    ylabel('mean signal correlation across songs')
    
    
    %using SignalCOrr and Distance, fit (exponential) relationship between them
    %to get Tau of spatial falloff of Signal Correlation
    figure;hold on
    line(distancebins_wide,Rmean_bydistance_wide,'color','m')
    line(distancebins_narrow,Rmean_bydistance_narrow,'color','b')
    axis tight
    set(gca,'YLim',[-0.2,1])
    title(BirdID)
    xlabel('distance between cluster')
    ylabel('mean signal correlation across songs')
    legend('WideSpikes','NarrowSpikes')
    
end

if do_fit ==1;
    maxD = 600; %don't try to fit too far out because much more noisy
    
    wide_fitinds = distancebins_wide<maxD;
    Dfit_wide = distancebins_wide(wide_fitinds);
    Rfit_wide = Rmean_bydistance_wide(wide_fitinds);
    fitobject_wide = fit(Dfit_wide,log(Rfit_wide)','poly1');
    SpaceConstant_wide = fitobject_wide.p1
    offset_wide = fitobject_wide.p2
    ci_wide = confint(fitobject_wide);
    ci_wide = ci_wide(:,1)
    
    narrow_fitinds = distancebins_narrow<maxD;
    Dfit_narrow = distancebins_narrow(narrow_fitinds);
    Rfit_narrow = Rmean_bydistance_narrow(narrow_fitinds);
    fitobject_narrow = fit(Dfit_narrow,log(Rfit_narrow)','poly1');
    SpaceConstant_narrow = fitobject_narrow.p1
    offset_narrow = fitobject_narrow.p2
    ci_narrow = confint(fitobject_narrow);
    ci_narrow = ci_narrow(:,1)
    
    
    if do_visualize ==1
        %check assumption of exponential
        figure;hold on
        scatter(distancebins_wide,log(Rmean_bydistance_wide),100,'b','fill')
        scatter(distancebins_narrow,log(Rmean_bydistance_narrow),100,'r','fill')
        title('assume exponential law (plot D by logR)')
        legend('Wide','Narrow')
        line(Dfit_wide,(fitobject_wide.p1*Dfit_wide)+fitobject_wide.p2,'color','k','LineWidth',3)
        line(Dfit_narrow,(fitobject_narrow.p1*Dfit_narrow)+fitobject_narrow.p2,'color','k','LineWidth',3)
        
        %check assumption of exponential
        figure;hold on
        scatter(Dfit_wide,log(Rfit_wide),100,'b','fill')
        scatter(Dfit_narrow,log(Rfit_narrow),100,'r','fill')
        title('assume exponential law (plot D by logR)')
        legend('Wide','Narrow')
        line(Dfit_wide,(fitobject_wide.p1*Dfit_wide)+fitobject_wide.p2,'color','k','LineWidth',3)
        line(Dfit_narrow,(fitobject_narrow.p1*Dfit_narrow)+fitobject_narrow.p2,'color','k','LineWidth',3)
        %
        % %check assumption of power law
        % figure;hold on
        % scatter(log(distancebins_wide),log(Rmean_bydistance_wide),100,'b','fill')
        % scatter(log(distancebins_narrow),log(Rmean_bydistance_narrow),100,'r','fill')
        % title('assume power law (plot logD by logR)')
        % legend('Wide','Narrow')
    end
end

%%
MeanRs_prestim = squeeze(nanmean(PreStimulusSignalCorr,3));

idx = true(size(MeanRs,1));
idx = ~(tril(idx));

R = MeanRs(idx);
D = PairWiseDistance(idx);

figure;scatter(D,R)


%% for each bird, plot the mean population response for each stimulus (wide versus narrow)

cd('/Users/kperks/mnt/cube/Ice/kperks/B1235/MetaMat/')

load('behavior.mat')
songinds = union(trainedSong_ind,novelSong_ind);


nstims = size(songinds,2);
matlist = {'TrialAveragedResponse_WideSpikes.mat','TrialAveragedResponse_NarrowSpikes.mat'};
popR_sum = [];
popR_mean = [];
for imat =1:2
    
    load(matlist{imat},'SignalFiltR','include_inds');
    Npop(imat) = size(include_inds,1);
%     if imat==1
%         Rmat(imat
%     end
%     if imat==2
%         
%     end
% end
% for imat = 1:2
    for istim = 1:nstims
        popR_sum(imat,songinds(istim),:) = squeeze(sum(SignalFiltR(:,songinds(istim),:),1));
        popR_mean(imat,songinds(istim),:) = squeeze(nanmean(SignalFiltR(:,songinds(istim),:),1));
    end
    
end

Ratio = squeeze(popR_mean(2,:,:)./popR_mean(1,:,:));
MagWide = squeeze(popR_mean(1,:,:));
MagNarrow = squeeze(popR_mean(2,:,:));
figure;
scatter(MagWide(:),Ratio(:))
title('dependence of I/E ratio on E')
xlabel('Magnitude mean Narrow Response')
ylabel('Narrow/Wide reponse ratio')
% figure;
% scatter(MagNarrow(:),Ratio(:))
MagBins = [0:0.001:max(MagWide(:))+0.01];
Rbinned = nan(1,size(MagBins,2)-1);
for ibin = 1:size(MagBins,2)-1
    maginds = intersect(find(MagWide>MagBins(ibin)),find(MagWide<MagBins(ibin+1)));
    Rbinned(ibin) = mean(Ratio(maginds));
    Mbinned(ibin) = mean(MagWide(maginds));
end
figure;
scatter(MagBins(2:end),Rbinned,400,'k','fill')
title('dependence of I/E ratio on E')
xlabel('Magnitude mean Narrow Response (mean by bin)')
ylabel('Narrow/Wide reponse ratio (mean by bin)')


figure;
hold on
for istim = 1:size(popR_sum,2)
    subplot(4,round(size(popR_sum,2)/4),istim);hold on
    plot(squeeze(popR_sum(1,istim,:))','color','b')
    plot(squeeze(popR_sum(2,istim,:))','color','r')
    axis tight
    ylabel(metatoes_songlist{songinds(istim)})
    
end

figure;
hold on
for istim = 1:size(popR_mean,2)
    subplot(4,round(size(popR_mean,2)/4),istim);hold on
    plot(squeeze(popR_mean(1,istim,:))','color','b')
    plot(squeeze(popR_mean(2,istim,:))','color','r')
    axis tight
    ylabel(metatoes_songlist{songinds(istim)})
    
end

idx = true(2);
idx = ~tril(idx);
popCorr = [];
popErr_norm = [];
popErr = [];
wideRnorm = [];
narrowRnorm = [];
Ratio = [];
MagWide = [];
for istim = 1:size(popR_mean,2)
    
   narrowR = squeeze(popR_mean(2,istim,:));
   narrowRnorm(istim,:) = narrowR ./ max(narrowR);
    
   wideR = squeeze(popR_mean(1,istim,:));
   wideRnorm(istim,:) = wideR./max(wideR);
   
    c = corrcoef(wideR,narrowR);
    popCorr(istim) = c(idx);
    
    popErr(istim,:) = narrowR - wideR;
    popErr_norm(istim,:) = narrowRnorm(istim,:) - wideRnorm(istim,:);
end

hfig = figure;
hold on
scatter(popCorr(novelSong_ind),popCorr(trainedSong_ind),400,'k','fill')
scatter(mean(popCorr(novelSong_ind)),mean(popCorr(trainedSong_ind)),400,'r')
xlabel('novel songs')
ylabel('trained songs')
title('correlation between mean population stim response')
SetAxisUnity(hfig)
lims = [0.5,1];
set(gca,'Xlim',lims,'YLim',lims)
% 
% figure;
% hold on
% plot(popErr_norm(novelSong_ind,:)','color','m')
% plot(popErr_norm(trainedSong_ind,:)','color','c')
% 
% novelErr = popErr_norm(novelSong_ind,:);
% N = novelErr(find(novelErr<0));
% trainedErr = popErr_norm(trainedSong_ind,:);
% T = trainedErr(find(trainedErr<0));
% round(sum(N))
% round(sum(T))
% 
% novelErr = popErr_norm(novelSong_ind,:);
% N = novelErr(find(novelErr>0));
% trainedErr = popErr_norm(trainedSong_ind,:);
% T = trainedErr(find(trainedErr>0));
% round(sum(N))
% round(sum(T))
% 

figure;
hold on
for istim = 1:size(popR_mean,2)
    subplot(4,round(size(popR_mean,2)/4),istim);hold on
    plot(narrowRnorm(istim,:)','color','r')
    plot(wideRnorm(istim,:)','color','b')
    axis tight
    ylabel(metatoes_songlist{songinds(istim)})
    
end


%% look at intracellular data in each bird

for icell = 1:size(IntracellularData,2)
    figure;hold on
    thiscell = IntracellularData{icell};
    for istim = 1:size(thiscell.stim_data,2)
        thisstim = thiscell.stim_data{istim};
        subplot(2,size(thiscell.stim_data,2)/2,istim)
        plot(mean(thisstim))
        axis tight
        ylabel(thiscell.stimnames{istim})
    end
end

%% plot extracellularly recorded population activity putatively around each of the intracellular sites
% cd('/Users/kperks/mnt/cube/Ice/kperks/B1235/MetaMat/')
intracellular_location = 1650;
breadth = 200;

load('behavior.mat')
songinds = union(trainedSong_ind,novelSong_ind);

% path_to_phy = 'Z:\Ice\kperks\B1112\klusta\phy041316\';
path_to_phy = '../klusta/phy041216/';
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
    load([path_to_phy sites{id} '/waveform_widths.mat'])
    thesesites = zeros(1,size(struct.cluster,2));
    thesesites(:) = id;
    siteID = [siteID, thesesites];
    clusterID = [clusterID, struct.cluster];
    isNarrow_TP = [isNarrow_TP, struct.narrow_trough_peak];
    isNarrow_HH = [isNarrow_HH, struct.narrow_half_height];

end
isNarrow = isNarrow_TP .* isNarrow_HH;

metatoes = load('metatoes.mat');
metatoesVars = metatoes;
vars = fieldnames(metatoes);
for ifield = 1:size(vars,1)
    s = [vars{ifield} '= metatoesVars.' vars{ifield} ';'];
    eval(s)
end

fs = metatoes{1}.fs;
stimendtime = metatoes{1}.stims{1}.stim_end_times-metatoes{1}.stims{1}.stim_start_times;
stimdur = stimendtime(1)/fs;

sigon = min(find(xtime_spikes>0));
sigoff = max(find(xtime_spikes<stimdur));

% calcstart = sigon + round((sigoff-sigon)/6);

% Good Quality
% From one penetration
qualinds = find(allsortclass == 2);
siteinds = find(allsiteID); %find(allsiteID == 1 | allsiteID == 2 | allsiteID == 3);

locationinds = intersect(find(allchanxz>(intracellular_location-breadth)),find(allchanxz<(intracellular_location+breadth)));

%%%%%%%%%%%%%%%%
% Wide Spike Clusters
%%%%%%%%%%%%%%%%
celltypeinds = find(isNarrow == 0);
include_inds = intersect(intersect(intersect(qualinds,siteinds),celltypeinds),locationinds);

nclusters = size(include_inds,1);

nstims = size(allSpikes_SynFilter,2);
ntrials = size(allSpikes_SynFilter,3);
nsamps = size(xtime_spikes,2);

TheseSpikes = squeeze(allSpikes_SynFilter(include_inds,:,:));
these_locations = allchanxz(include_inds);

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

%%%%%%%%%%%%%%%%
% Narrow Spike Clusters
%%%%%%%%%%%%%%%%
celltypeinds = find(isNarrow == 1);
include_inds = intersect(intersect(intersect(qualinds,siteinds),celltypeinds),locationinds);

nclusters = size(include_inds,1);

TheseSpikes = squeeze(allSpikes_SynFilter(include_inds,:,:));
these_locations = allchanxz(include_inds);

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
            SignalFiltR_narrow(isite,istim,:) = mean(TrialFiltR(:,sigon:sigoff));
        end
end

figure;hold on
for istim = 1:nstims
    subplot(4,round(nstims/4),istim);hold on
    line(xtime_spikes(sigon:sigoff),squeeze(mean(SignalFiltR_wide(:,istim,:),1)),'color','b')
    line(xtime_spikes(sigon:sigoff),squeeze(mean(SignalFiltR_narrow(:,istim,:),1)),'color','r')
    ylabel(metatoes_songlist{istim})
    axis tight
end


figure;hold on
for istim = 1:nstims
    subplot(4,round(nstims/4),istim);hold on
    wide_mean = squeeze(mean(SignalFiltR_wide(:,istim,:),1));
    line(xtime_spikes(sigon:sigoff),wide_mean ./ max(wide_mean),'color','b')
    narrow_mean = squeeze(mean(SignalFiltR_narrow(:,istim,:),1));
    line(xtime_spikes(sigon:sigoff),narrow_mean./max(narrow_mean),'color','r')
    ylabel(metatoes_songlist{istim})
    axis tight
end


