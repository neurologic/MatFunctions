load('TrialAveragedResponse_WideSpikes.mat','SignalFiltR')
load('SpatialSignalCorr_WideSpikes.mat', 'PairWiseDistance','SignalCorr','include_inds');

nclusters = size(include_inds,1);

nstims = size(SignalFiltR,2);
% ntrials = size(allSpikes_SynFilter,3);
nsamps = size(SignalFiltR,3);

GeoMean = zeros(nclusters,nclusters,nstims);
GeoMeanRate = zeros(nclusters,nclusters,nstims);
for isite1 = 1:nclusters
    siteprogress_countdown = nclusters - isite1
    
    %     fileID = fopen('matlog.txt','w');
    %     t = datestr(datetime('now'));
    %     fprintf(fileID,'%s %i %s\n','Wide_countdown',siteprogress_countdown,t);
    %     fclose(fileID);
    
    for isite2 = isite1:nclusters
        
        for istim = 1:nstims
            R1 = squeeze(SignalFiltR(isite1,istim,:));
            R2 = squeeze(SignalFiltR(isite2,istim,:));
            
            
%             GeoMean(isite1,isite2,istim) = sum(geomean([R1';R2'])');
            GeoMeanRate(isite1,isite2,istim) = (geomean([sum(R1),sum(R2)]));
            
        end
    end
end

% save('SpatialSignalCorr_WideSpikes.mat','PairWiseDistance','SignalCorr','include_inds','GeoMean')
load('behavior.mat')
songinds = union(trainedSong_ind,novelSong_ind);

MeanGeoMean = squeeze(mean(GeoMean(:,:,songinds),3));
MeanRS = squeeze(mean(SignalCorr(:,:,songinds),3));
%get indices of upper triangle of corr matrix
idx = true(size(MeanGeoMean,1));
idx = ~(tril(idx));
% get vector for upper triangle of SignalCorr and Distance
M = MeanGeoMean(idx);
R = MeanRS(idx);

figure;scatter(M,abs(R))
ylabel('abs(mean pairwise Signal Correlation across stimuli)')
xlabel('mean pairwise Geometric mean between responses')

figure;hold on
ncol = 4
for istim = 1:nstims
    subplot(round(nstims/ncol),ncol,istim)
    stimR = squeeze(SignalCorr(:,:,istim));
    stimG = squeeze(GeoMeanRate(:,:,istim));
    scatter(stimR(idx),stimG(idx))
end


