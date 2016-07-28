cd('/mnt/cube/Ice/kperks/B1112/MetaMat/')
currentpath = pwd
% addpath('/mnt/cube/Ice/kperks/MatFunctions/');
%%%%%%%%%%%%%%%%
%load and set up data for analysis
%%%%%%%%%%%%%%%%


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
siteinds = find(allsiteID == 1 | allsiteID == 2 | allsiteID == 3);find(allsiteID); %

%%%%%%%%%%%%%%%%
% All Spike Clusters
%%%%%%%%%%%%%%%%
% celltypeinds = find(isNarrow == 0);
% include_inds = intersect(intersect(qualinds,siteinds),celltypeinds);
include_inds = intersect(qualinds,siteinds);

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
    fprintf(fileID,'%s %i %s\n','Prep countdown',siteprogress_countdown,t);
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
save('TrialAveragedResponse_AllSpikes_Sites123.mat','SignalFiltR','include_inds')

SignalCorr = nan(nclusters,nclusters,nstims);
PairWiseDistance = nan(nclusters,nclusters);

%need specific index for corr matrix because sometimes all NaN
idx = true(2);
idx = ~(tril(idx));

for isite1 = 1:nclusters
    siteprogress_countdown = nclusters - isite1;

    fileID = fopen('matlog.txt','w');
    t = datestr(datetime('now'));
    fprintf(fileID,'%s %i %s\n','Corr calc countdown',siteprogress_countdown,t);
    fclose(fileID);

    for isite2 = isite1:nclusters
        
        for istim = 1:nstims
          

            c = corrcoef(squeeze(SignalFiltR(isite1,istim,:)),squeeze(SignalFiltR(isite2,istim,:)));
           
            SignalCorr(isite1,isite2,istim) = c(idx);
            
            
        end
        PairWiseDistance(isite1,isite2) = abs(these_locations(isite2) - these_locations(isite1));
    end
end

locations = allchanxz(include_inds);

[SoretedLocations,SortLocationInds] = sort(locations);
SortedRhoSignal = SignalCorr(SortLocationInds,SortLocationInds,:);
% 
% fullmat = [];
% fullmat(1,:,:) = mean(SignalCorr,3);
% fullmat(2,:,:) = mean(SignalCorr,3)';
% symmetricmat = squeeze(nanmean(fullmat,1));
% symmetricmat(eye(size(symmetricmat))~=0) = nan;
% symmetricmat(isnan(symmetricmat)) = 0;
% figure;
% imagesc(symmetricmat(SortLocationInds,SortLocationInds));
% n = 20;              %// number of colors
% R = linspace(0,1,n);  %// Red from 1 to 0
% G = linspace(0,1,n);   %// Green 1 to 0
% B = linspace(1,0,n);  %// Blue from 0 to 1
% colormap( [R(:), G(:), B(:)] );  %// create colormap
% caxis([-1, 1])
% colorbar
% xtickinds = get(gca,'XTickLabel');
% xtickinds = cellfun(@str2num, xtickinds);
% ytickinds = get(gca,'YTickLabel');
% ytickinds = cellfun(@str2num, ytickinds);
% set(gca,'XTick',xtickinds,'XTickLabel',SoretedLocations(xtickinds),'YTick',ytickinds,'YTickLabel',SoretedLocations(ytickinds))
% 
% xlims = get(gca,'XLim');
% ylims = get(gca,'YLim');

save('SpatialResponseSimilarity_AllSpikes_Sites123.mat','SortedRhoSignal','SoretedLocations','SortLocationInds','PairWiseDistance','SignalCorr','include_inds')



