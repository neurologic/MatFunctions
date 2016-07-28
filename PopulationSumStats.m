
% cd('/Users/kperks/mnt/cube/Ice/kperks/B1112/MetaMat')
cd('Z:\Ice\kperks\B1112\MetaMat')
%%

metatoes = load('metatoes.mat');
metatoesVars = metatoes;
vars = fieldnames(metatoes);
for ifield = 1:size(vars,1)
    s = [vars{ifield} '= metatoesVars.' vars{ifield} ';'];
    eval(s)
end

sigon = min(find(xtime_spikes>0));
sigoff = max(find(xtime_spikes<6));
%%

siteID_list = unique(allsiteID);

do_trialshuffle_list = repmat([0,1],1,size(siteID_list,2));
shufflename = {'Normal','Shuffle'};

datamat_list = []

restrict_clusterQ = 2;
%%
for isite = 1:size(siteID_list,2)
    for idoshuffle = 1:2
        do_trialshuffle = do_trialshuffle_list(((isite*2)-1) + (idoshuffle -1));
        
        datamat_list_name = ['PopSumData_' shufflename{idoshuffle} '_site' num2str(siteID_list(isite)) '.mat'];
        datamat_list{(isite+do_trialshuffle)} = datamat_list_name;

        
        %defaults inds
        sortclassinds = [1:size(allsortclass,2)];
        siteinds = [1:size(allsiteID,2)];
        
        sortclassinds = find(allsortclass == restrict_clusterQ);
        siteinds = find(allsiteID == siteID_list(isite));
        include_inds = intersect(sortclassinds,siteinds);
        
        nclusters = size(include_inds,1);
        nstims = size(allSpikes_SynFilter,2);
        ntrials = size(allSpikes_SynFilter,3);
        nsamps = size(xtime_spikes,2);
        
        ResponseSumFiltR = zeros(nstims,ntrials,nsamps);
        
        for istim = 1:size(metatoes{1}.stims,1)
            
            TrialFiltR = zeros(nclusters,nsamps);
            
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
                
            end
            
            ResponseSumFiltR(istim,:,:) = meanpsth;
            
        end
        save(datamat_list_name, 'ResponseSumFiltR')
        
    end
end

save('PopSum_datamat_list.mat','datamat_list')

%%
datamat_list = [];
for isite = 1:size(siteID_list,2)
    for idoshuffle = 1:2
        do_trialshuffle = do_trialshuffle_list(isite + idoshuffle -1);
        
        datamat_list_name = ['PopSumData_' shufflename{idoshuffle} '_site' num2str(siteID_list(isite)) '.mat'];
        datamat_list{(isite*2)-1+do_trialshuffle} = datamat_list_name;
        
    end
end
save('PopSum_datamat_list.mat','datamat_list')
%% Stats of population sum activity normal versus trial-shuffled data per cluster

% load('PopSum_datamat_list.mat','datamat_list')
%vars in datamat_list to analyze
% datamat_list = {'DataNorm_site1','DataShuffle_site1','DataNorm_site2','DataShuffle_site2'}

DataTable = cell(2,2,6,size(datamat_list,2)/2,nstims);

for isite = 1:size(datamat_list,2);
    
    %shuffle or normal data (in data table, normal data is row 1 and shuffled data is row 2)
    if ~isempty(regexp(datamat_list{isite},'Normal'));
        DataRowInd = 1;
    end
    if ~isempty(regexp(datamat_list{isite},'Shuffle'));
        DataRowInd = 2;
    end
    
    a = regexp(datamat_list{isite},'_');
    b = regexp(datamat_list{isite},'.mat');
    siteID = str2num(datamat_list{isite}(a(2)+5:b-1));
    
    load(datamat_list{isite},'ResponseSumFiltR');
    
    for istim = 5:nstims% 1:4%nstims
        
        %get trials in response to one pseudo song
        stimdata = squeeze(ResponseSumFiltR(istim,:,:));
        
        %normalize to max across all trials
        stimdata = stimdata ./ max(max(stimdata));
        
        startind = 225;
        PreStimData = stimdata(:,startind:sigon);
        StimData = stimdata(:,sigon:sigoff);
        ColumnData_List = {'PreStimData','StimData'};
        
        % {x,1} PreStim ;x = DataRowInd
        % {x,2} Stim ;x = DataRowInd
        
        for DataColumnInd = 1:2
            
            s = ['stimdata = ' ColumnData_List{DataColumnInd} ';'];
            eval(s)
            
            %calculate the signal mean
            trialavg = mean(stimdata,1);
            
            %calcualte the residuals
            trialresid = stimdata - repmat(trialavg,size(stimdata,1),1);
            
            % 1:skew Vm all
            DataTable{DataRowInd,DataColumnInd,1,siteID,istim} = skewness(stimdata(:));
            
            % 2:skew Vm residuals
            DataTable{DataRowInd,DataColumnInd,2,siteID,istim} = skewness(trialresid(:));
            
            % 3:skew Vm signal
            DataTable{DataRowInd,DataColumnInd,3,siteID,istim} = skewness(trialavg(:));
            
            % 4:trial var Vm
            DataTable{DataRowInd,DataColumnInd,4,siteID,istim} = var(stimdata,1);
            
            % 5:trial-avg Vm
            DataTable{DataRowInd,DataColumnInd,5,siteID,istim} = trialavg;
            
            % 6:all trials Vm
            DataTable{DataRowInd,DataColumnInd,5,siteID,istim} = stimdata;
        end
        
    end
    
end

%%
c = colormap(parula(6));

hfig = figure,hold on; 
for isite = 1:size(DataTable,4)
prestim = cell2mat(squeeze(DataTable(1,1,2,isite,:)));
stim = cell2mat(squeeze(DataTable(2,1,2,isite,:)));
scatter (prestim,stim,100,c(isite,:),'fill')
end
SetAxisUnity(hfig)
title('prestim Vm residual skew')
xlabel('normal')
ylabel('shuffled')

hfig = figure,hold on; 
for isite = 1:size(DataTable,4)
prestim = cell2mat(squeeze(DataTable(1,2,2,isite,:)));
stim = cell2mat(squeeze(DataTable(2,2,2,isite,:)));
scatter (prestim(:),stim(:),100,c(isite,:),'fill')
end
SetAxisUnity(hfig)
title('stim Vm residual skew')
xlabel('normal')
ylabel('shuffled')


normal = mean((cell2mat(squeeze(DataTable(1,2,2,:,:))) ./ cell2mat(squeeze(DataTable(1,1,2,:,:)))),2);
shuffle = mean((cell2mat(squeeze(DataTable(2,2,2,:,:))) ./ cell2mat(squeeze(DataTable(2,1,2,:,:)))),2);
hfig = figure;
scatter(normal,shuffle,100,'k','fill')
title('Vm residual skew')
xlabel('stim/prestim under normal condition')
ylabel('stim/prestim under trial shuffle condition')
SetAxisUnity(hfig)

%%
c = colormap(parula(6));

hfig = figure,hold on; 
for isite = 1:size(DataTable,4)
prestim = cell2mat(squeeze(DataTable(1,1,1,isite,:)));
stim = cell2mat(squeeze(DataTable(2,1,1,isite,:)));
scatter (prestim,stim,100,c(isite,:),'fill')
end
SetAxisUnity(hfig)
title('prestim Vm(total) skew')
xlabel('normal')
ylabel('shuffled')

hfig = figure,hold on; 
for isite = 1:size(DataTable,4)
prestim = cell2mat(squeeze(DataTable(1,2,1,isite,:)));
stim = cell2mat(squeeze(DataTable(2,2,1,isite,:)));
scatter (prestim(:),stim(:),100,c(isite,:),'fill')
end
SetAxisUnity(hfig)
title('stim Vm(total) skew')
xlabel('normal')
ylabel('shuffled')


normal = mean((cell2mat(squeeze(DataTable(1,2,1,:,:))) ./ cell2mat(squeeze(DataTable(1,1,1,:,:)))),2);
shuffle = mean((cell2mat(squeeze(DataTable(2,2,1,:,:))) ./ cell2mat(squeeze(DataTable(2,1,1,:,:)))),2);
hfig = figure;
scatter(normal,shuffle,100,'k','fill')
title('Vm(total) skew')
xlabel('stim/prestim under normal condition')
ylabel('stim/prestim under trial shuffle condition')
SetAxisUnity(hfig)

%%
c = colormap(parula(6));

hfig = figure,hold on; 
for isite = 1:size(DataTable,4)
prestim = mean(cell2mat(squeeze(DataTable(1,1,4,isite,:))),2);
stim = mean(cell2mat(squeeze(DataTable(2,1,4,isite,:))),2);
scatter (prestim,stim,100,c(isite,:),'fill')
end
SetAxisUnity(hfig)
title('pre-stim trial var Vm')
xlabel('normal')
ylabel('shuffled')

hfig = figure,hold on; 
for isite = 1:size(DataTable,4)
prestim = mean(cell2mat(squeeze(DataTable(1,2,4,isite,:))),2);
stim = mean(cell2mat(squeeze(DataTable(2,2,4,isite,:))),2);
scatter (prestim(:),stim(:),100,c(isite,:),'fill')
end
SetAxisUnity(hfig)
title('stim trial var Vm')
xlabel('normal')
ylabel('shuffled')


normal = mean(mean(cell2mat(squeeze(DataTable(1,2,4,:,:))),2) ./ mean(cell2mat(squeeze(DataTable(1,1,4,:,:))),2),2);
shuffle = mean(mean(cell2mat(squeeze(DataTable(2,2,4,:,:))),2) ./ mean(cell2mat(squeeze(DataTable(2,1,4,:,:))),2),2);
hfig = figure;
scatter(normal,shuffle,100,'k','fill')
title('trial var Vm')
xlabel('stim/prestim under normal condition')
ylabel('stim/prestim under trial shuffle condition')
SetAxisUnity(hfig)

%% Signal Correlation Per Cluster against distance between clusters
group1 = find(allsiteID == 1 | allsiteID == 2 | allsiteID ==3);
siteinds = group1;

restrict_clusterQ = 2;
sortclassinds = find(allsortclass == restrict_clusterQ);

include_inds = intersect(sortclassinds,siteinds);
nclusters = size(include_inds,1);

nstims = size(allSpikes_SynFilter,2);
ntrials = size(allSpikes_SynFilter,3);
nsamps = size(xtime_spikes,2);

TheseSpikes = squeeze(allSpikes_SynFilter(include_inds,:,:));

these_locations = allchanxz(include_inds);

SignalCorr = nan(nclusters,nclusters,nstims);
PairWiseDistance = nan(nclusters,nclusters);
for isite1 = 1:nclusters
    for isite2 = isite1:nclusters
        

        for istim = 1:size(metatoes{1}.stims,1)
            
            TrialFiltR = zeros(2,ntrials,nsamps);
            for itrial = 1:ntrials
                trialfun1 = TheseSpikes{isite1,istim,itrial};
                if ~isempty(trialfun1)
                    TrialFiltR(1,itrial,:) = trialfun1(xtime_spikes); %for this cluster on this trial
                end
                 trialfun2 = TheseSpikes{isite2,istim,itrial};
                if ~isempty(trialfun2)
                    TrialFiltR(2,itrial,:) = trialfun2(xtime_spikes); %for this cluster on this trial
                end
                
            end
            c = corrcoef(mean(squeeze(TrialFiltR(1,:,sigon:sigoff))),mean(squeeze(TrialFiltR(2,:,sigon:sigoff))));
            SignalCorr(isite1,isite2,istim) = c(~triu(c));
            
            PairWiseDistance(isite1,isite2) = abs(these_locations(isite2) - these_locations(isite1));
        end
    end
end

   


%%

load('ExtracellularSignalCorrAnalysis.mat','PairWiseDistance', 'SignalCorr', 'PreStimulusSignalCorr')

MeanRs = squeeze(nanmean(SignalCorr,3));

idx = true(size(MeanRs,1)); 
idx = ~(tril(idx));

R = MeanRs(idx);
D = PairWiseDistance(idx);

distancebins = unique(D);

cols = colormap(hsv(size(distancebins,1)));

figure;hold on
for iD = 1:size(distancebins,1)

    theseinds = find(D == distancebins(iD));
    
    scatter(D(theseinds),R(theseinds),10,cols(iD,:))
    scatter(distancebins(iD),nanmean(R(theseinds)),200,'k','fill')
    
    Rmean_bydistance(iD) = nanmean(R(theseinds));

end

figure;hold on
line(distancebins,Rmean_bydistance,'color','m')


MeanRs_prestim = squeeze(nanmean(PreStimulusSignalCorr,3));

idx = true(size(MeanRs,1)); 
idx = ~(tril(idx));

R = MeanRs(idx);
D = PairWiseDistance(idx);

figure;scatter(D,R)

%% get mean psth per electrode port and then get correlation across ports and plot by distance
load('TrialAveraged_Sites123_ClusterQuality2.mat','TrialAveraged_Stimulus','PairWiseDistance')
load('ExtracellularSignalCorrAnalysis.mat','PairWiseDistance')

locationbins = unique(these_locations);

PortSum = zeros(size(locationbins,1),nstims,size(TrialAveraged_Stimulus,3));


for iport = 1:size(locationbins,1)
    thisportinds = find(these_locations == locationbins(iport));
    
    PortSum(iport,:,:) = squeeze(mean(TrialAveraged_Stimulus(thisportinds,:,:),1));

end

idx = true(2); 
idx = ~(tril(idx));
NetPortCorr = nan(size(PortSum,1),size(PortSum,1),nstims);
PortDistances = nan(size(PortSum,1),size(PortSum,1));
for iport1 = 1:size(PortSum,1)
    for iport2 = iport1:size(PortSum,1)
        
        PortDistances(iport1,iport2) = abs(locationbins(iport1)-locationbins(iport2));

        for istim = 1:size(PortSum,2)
            c = corrcoef(PortSum(iport1,istim,:),PortSum(iport2,istim,:));
            NetPortCorr(iport1,iport2,istim) = c(idx);
        end
    end 
end


MeanRs_port = squeeze(nanmean(NetPortCorr,3));

idx = true(size(MeanRs_port,1)); 
idx = ~(tril(idx));

R = MeanRs_port(idx);
D = PortDistances(idx);

distancebins = unique(D);

cols = colormap(hsv(size(distancebins,1)));

figure;hold on
for iD = 1:size(distancebins,1)

    theseinds = find(D == distancebins(iD));
    
%     scatter(D(theseinds),R(theseinds),10,cols(iD,:))
%     scatter(distancebins(iD),nanmean(R(theseinds)),200,'k','fill')
    
    Rmean_bydistance_port(iD) = nanmean(R(theseinds));
end

% figure;hold on
line(distancebins,Rmean_bydistance_port,'color','b')

%%
load('ExtracellularResidualCorrAnalysis_site1.mat', 'ResidualsCorr','PreStimResidualsCorr')
load('TrialResiduals_Site1_ClusterQuality2.mat','PairWiseDistance')

MeanRs = squeeze(nanmean(ResidualsCorr,4));
MeanRs = squeeze(nanmean(MeanRs,3));

idx = true(size(MeanRs,1)); 
idx = ~(tril(idx));

R = MeanRs(idx);
D = PairWiseDistance(idx);

distancebins = unique(D);

cols = colormap(hsv(size(distancebins,1)));

Rmean_bydistance = nan(1,size(distancebins,1));
figure;hold on
for iD = 1:size(distancebins,1)

    theseinds = find(D == distancebins(iD));
    
    scatter(D(theseinds),R(theseinds),10,cols(iD,:))
    scatter(distancebins(iD),nanmean(R(theseinds)),200,'k','fill')
    
    Rmean_bydistance(iD) = nanmean(R(theseinds));

end

figure;hold on
line(distancebins,Rmean_bydistance,'color','m')

%% need to get residuals correlation by bin
%%%%% try both time series corr within each bin and try 

%%%%%%%%%% wait... is correlating time series of risiduals even the right
%%%%%%%%%% way to ask? 

% normally when look at noise corr one would look at correlating trial
% variance for a given sample/stimulus/bin

%..... may just be two different questions... both valid to ask?
% load('TrialResiduals_Site1_ClusterQuality2.mat','TrialResiduals_Stimulus')


load('SignalTrialCorrs_bin100ms_site1_ClusterQuality2.mat','SignalTrialCorr')
load('SignalGeoMean_bin100ms_site1_ClusterQuality2.mat','SignalGeoMean')

% FOR EACH cell pair... get correlation between trial corr and geo mean
% across signals
idc = true(2);
idc = ~tril(idc);
nclusters = size(SignalTrialCorr,1);
nstims = size(SignalTrialCorr,3);
nbins = size(SignalTrialCorr,4);

TrialCorr_V_GeoMean = nan(nclusters,nclusters);
idx = true(nclusters);
idx = ~tril(idx);

figure
hold on
for isite1 = 1:nclusters
   for isite2 = (isite1+1):nclusters 
%        isite1 = 20; isite2 = 50;
       thispair_corrs = squeeze(SignalTrialCorr(isite1,isite2,:,:));
       thispair_corrs = thispair_corrs(:);
       thispair_geomean = squeeze(SignalGeoMean(isite1,isite2,:,:));
       thispair_geomean = thispair_geomean(:) / max(max(thispair_geomean));
       
       %sort geomean into bins and get mean geomean and mean this_pair
       %corrs for corresponding stimuli
       
       bins = [0:0.1:1];
       for ibin = 1:size(bins,2)-1
          indslow = find(thispair_geomean>bins(ibin));
          indshigh = find(thispair_geomean<bins(ibin+1));
          inds = intersect(indslow,indshigh);
          binned_geomean(isite1,isite2,ibin) = nanmean(thispair_geomean(inds));
          binned_corr(isite1,isite2,ibin) = nanmean(thispair_corrs(inds));
       end
       
       
       scatter(binned_geomean(isite1,isite2,:),binned_corr(isite1,isite2,:))
       c = corrcoef(binned_corr(isite1,isite2,:),binned_geomean(isite1,isite2,:));
       TrialCorr_V_GeoMean(isite1,isite2) = c(idc);
   end
end

n = histc(TrialCorr_V_GeoMean(:),[-1:0.1:1]);
n = n/size(TrialCorr_V_GeoMean(:),1);
figure;stairs([-1:0.1:1],n)

%now bin by trial correlation and don't normalize geometric mean
TrialCorr_V_GeoMean = nan(nclusters,nclusters);
idx = true(nclusters);
idx = ~tril(idx);
figure
hold on
for isite1 = 1:nclusters
   for isite2 = (isite1+1):nclusters 
%        isite1 = 20; isite2 = 50;
       thispair_corrs = squeeze(SignalTrialCorr(isite1,isite2,:,:));
       thispair_corrs = thispair_corrs(:);
       thispair_geomean = squeeze(SignalGeoMean(isite1,isite2,:,:));
       thispair_geomean = thispair_geomean(:);% / max(max(thispair_geomean));
       
       %sort geomean into bins and get mean geomean and mean this_pair
       %corrs for corresponding stimuli
       
       bins = [-1:0.2:1];
       for ibin = 1:size(bins,2)-1
          indslow = find(thispair_corrs>bins(ibin));
          indshigh = find(thispair_corrs<bins(ibin+1));
          inds = intersect(indslow,indshigh);
          binned_geomean(isite1,isite2,ibin) = nanmean(thispair_geomean(inds));
          binned_corr(isite1,isite2,ibin) = nanmean(thispair_corrs(inds));
       end
       
       
       scatter(binned_geomean(isite1,isite2,:),binned_corr(isite1,isite2,:))
       c = corrcoef(binned_corr(isite1,isite2,:),binned_geomean(isite1,isite2,:));
       TrialCorr_V_GeoMean(isite1,isite2) = c(idc);
   end
end
n = histc(TrialCorr_V_GeoMean(:),[-1:0.1:1]);
n = n/size(TrialCorr_V_GeoMean(:),1);
figure;stairs([-1:0.1:1],n)


%%
load('TrialResiduals_Site1_ClusterQuality2.mat','PairWiseDistance')

figure;hold on
for isite1 = 1:nclusters
    for isite2 = (isite1+1):nclusters
        thispair_corrs = squeeze(SignalTrialCorr(isite1,isite2,:,:));
        thispair_corrs = thispair_corrs(:);
        
        scatter(max(thispair_corrs), PairWiseDistance(isite1,isite2));
        
    end
end

MeanRs = squeeze(nanmean(SignalTrialCorr,4));
MeanRs = squeeze(nanmean(MeanRs,3));

idx = true(size(MeanRs,1)); 
idx = ~(tril(idx));

R = MeanRs(idx);
D = PairWiseDistance(idx);

distancebins = unique(D);

cols = colormap(hsv(size(distancebins,1)));

Rmean_bydistance = nan(1,size(distancebins,1));
figure;hold on
for iD = 1:size(distancebins,1)

    theseinds = find(D == distancebins(iD));
    
    scatter(D(theseinds),R(theseinds),10,cols(iD,:))
    scatter(distancebins(iD),nanmean(R(theseinds)),200,'k','fill')
    
    Rmean_bydistance(iD) = nanmean(R(theseinds));

end


       

%% for a single site, trial correlations (correlation in residuals on each trial) for each Signal 
group1 = find(allsiteID == 2);% | allsiteID == 2 | allsiteID ==3);
siteinds = group1;

restrict_clusterQ = 2;
sortclassinds = find(allsortclass == restrict_clusterQ);

include_inds = intersect(sortclassinds,siteinds);
nclusters = size(include_inds,1);

nstims = size(allSpikes_SynFilter,2);
ntrials = size(allSpikes_SynFilter,3);
nsamps = size(xtime_spikes,2);

TheseSpikes = squeeze(allSpikes_SynFilter(include_inds,:,:));

these_locations = allchanxz(include_inds);

MeanTrialCorr = nan(nclusters,nclusters,nstims);
PairWiseDistance = nan(nclusters,nclusters);
for isite1 = 1:nclusters
    for isite2 = isite1:nclusters
        

        for istim = 1:size(metatoes{1}.stims,1)
            
            TrialFiltR = zeros(2,ntrials,nsamps);
            for itrial = 1:ntrials
                trialfun1 = TheseSpikes{isite1,istim,itrial};
                if ~isempty(trialfun1)
                    TrialFiltR(1,itrial,:) = trialfun1(xtime_spikes); %for this cluster on this trial
                end
                 trialfun2 = TheseSpikes{isite2,istim,itrial};
                if ~isempty(trialfun2)
                    TrialFiltR(2,itrial,:) = trialfun2(xtime_spikes); %for this cluster on this trial
                end
                
            end
            c = corrcoef(mean(squeeze(TrialFiltR(1,:,sigon:sigoff))),mean(squeeze(TrialFiltR(2,:,sigon:sigoff))));
            
            
            MeanTrialCorr(isite1,isite2,istim) = c(~triu(c));
            
            PairWiseDistance(isite1,isite2) = abs(these_locations(isite2) - these_locations(isite1));
        end
    end
end

%%

% for a subset of clsuters, 
% eliminate first 1/6 of response (1sec in good expts, 1.x in expts with
% mis-sampled rate
%for one subset of stimuli...

% calc mean and variance per stimulus
sigon = min(find(xtime_spikes>0));
sigoff = max(find(xtime_spikes<6));

calcstart = sigon + round((sigoff-sigon)/6); %size(TrialAvg_Stimulus,3)/6; %1/6 of way into 6-motif pseudo song

stimsize = floor(0.04*11025); %40msec

bins = [calcstart:stimsize:sigoff];

load('PopSumData_Normal_site1.mat')
ResponseSumFiltR_Normal = ResponseSumFiltR;

%novel songs
data = [];
stimind = 1;
for isong = 1:4
    thissong = squeeze(ResponseSumFiltR_Normal(isong,:,:));
    for ibin = 1:size(bins,2)-1
        data(stimind,:,:) = thissong(:,bins(ibin):(bins(ibin+1)-1));
        stimind = stimind+1;
    end
end

[m, v] = KLprep_Univariate(data);
kl_pairwise_NN = KLpairwise_Univariate(m , v);

%trained songs
data = [];
stimind = 1;
for isong = 5:8
    thissong = squeeze(ResponseSumFiltR_Normal(isong,:,:));
    for ibin = 1:size(bins,2)-1
        data(stimind,:,:) = thissong(:,bins(ibin):(bins(ibin+1)-1));
        stimind = stimind+1;
    end
end

[m, v] = KLprep_Univariate(data);
kl_pairwise_TN = KLpairwise_Univariate(m , v);
% figure;imagesc(KL_NetPseudoSyn_Trained_TrialNormal)
% title('KL_NetPseudoSyn_Trained_TrialNormal','Interpreter','none')


load('PopSumData_Shuffle_site1.mat')
ResponseSumFiltR_Shuffle = ResponseSumFiltR;

%novel songs
data = [];
stimind = 1;
for isong = 1:4
    thissong = squeeze(ResponseSumFiltR_Shuffle(isong,:,:));
    for ibin = 1:size(bins,2)-1
        data(stimind,:,:) = thissong(:,bins(ibin):(bins(ibin+1)-1));
        stimind = stimind+1;
    end
end

[m, v] = KLprep_Univariate(data);
kl_pairwise_NS = KLpairwise_Univariate(m , v);

%trained songs
data = [];
stimind = 1;
for isong = 5:8
    thissong = squeeze(ResponseSumFiltR_Shuffle(isong,:,:));
    for ibin = 1:size(bins,2)-1
        data(stimind,:,:) = thissong(:,bins(ibin):(bins(ibin+1)-1));
        stimind = stimind+1;
    end
end

[m, v] = KLprep_Univariate(data);
kl_pairwise_TS = KLpairwise_Univariate(m , v);
% figure;imagesc(KL_NetPseudoSyn_Trained_TrialNormal)
% title('KL_NetPseudoSyn_Trained_TrialNormal','Interpreter','none')

%%

nsubsamps = 1000;
nreps = 1000;

NN_sumKL = KLsum_Univariate(kl_pairwise_NN,nsubsamps,nreps);
TN_sumKL = KLsum_Univariate(kl_pairwise_TN,nsubsamps,nreps);
NS_sumKL = KLsum_Univariate(kl_pairwise_NS,nsubsamps,nreps);
TS_sumKL = KLsum_Univariate(kl_pairwise_TS,nsubsamps,nreps);

figure;hold on
[h,stats_TN] = cdfplot(TN_sumKL);
[h,stats_TS] = cdfplot(TS_sumKL);
[h,stats_NN] = cdfplot(NN_sumKL);
[h,stats_NS] = cdfplot(NS_sumKL);
legend('Trained TrialNormal','Trained TrialShuffle','Novel TrialNormal','Novel TrialShuffle')
xlabel(['sum KLdivergence for subset ' num2str(nsubsamps) ' stims; run ' num2str(nreps) ' times'])

figure;hold on
[h,stats_TN] = cdfplot(TN_sumKL);
[h,stats_NN] = cdfplot(NN_sumKL);
legend('trained','novel')

TrialShuffleChangeTrained = (TS_sumKL ./ TN_sumKL);
TrialShuffleChangeNovel = (NS_sumKL ./ NN_sumKL) ;
figure;hold on
cdfplot(TrialShuffleChangeNovel)
cdfplot(TrialShuffleChangeTrained)
legend('KLsum trial shuffle change Novel', 'KLsum trial shuffle change Trained')

edges = [0:0.05:4];
n_Novel = histc(TrialShuffleChangeNovel,edges)./size(TrialShuffleChangeNovel,2);
n_Trained = histc(TrialShuffleChangeTrained,edges)./size(TrialShuffleChangeTrained,2);
figure;hold on
stairs(edges,n_Novel,'color','k','LineWidth',5)
stairs(edges,n_Trained,'color','b','LineWidth',3)
legend('KLsum trial shuffle change Novel', 'KLsum trial shuffle change Trained')
h_TrailCorrEffect_NvsT = kstest2(TrialShuffleChangeNovel,TrialShuffleChangeTrained)

%%
idx = true(nstimbins,nstimbins);
idx = ~tril(idx);
edgesmax = round(max([max(TN_pairs),max(NN_pairs),max(NS_pairs),max(TS_pairs)]));
edgesmin = 0;
edges = [edgesmin:0.5:30];%edgesmax+mod(edgesmax,round(edgesmax/1000))];
n_TN = histc(TN_pairs,edges);
n_NN = histc(NN_pairs,edges);
n_TS = histc(TS_pairs,edges);
n_NS = histc(NS_pairs,edges);
figure;hold on
line(edges,n_TN,'color',[1,0,0],'LineWidth',4)
line(edges,n_TS,'color',[0.8,0,0],'LineWidth',4)
line(edges,n_NN','color',[0,0,1],'LineWidth',4)
line(edges,n_NS,'color',[0,0,0.8],'LineWidth',4)
legend('Trained TrialNormal','Trained_TrialShuffle','Novel_TrialNormal','Novel_TrialShuffle')

[n_TN,centersTN] = hist(TN_pairs);
[n_NN,centersNN] = hist(NN_pairs);
[n_TS,centersTS] = hist(TS_pairs);
[n_NS,centersNS] = hist(NS_pairs);
figure;hold on
stairs(centersTN,n_TN,'color',[1,0,0],'LineWidth',4)
stairs(centersTS,n_TS,'color',[0.8,0,0],'LineWidth',4)
stairs(centersNN,n_NN','color',[0,0,1],'LineWidth',4)
stairs(centersNS,n_NS,'color',[0,0,0.8],'LineWidth',4)
legend('Trained TrialNormal','Trained_TrialShuffle','Novel_TrialNormal','Novel_TrialShuffle')

figure;hold on
[h,stats_TN] = cdfplot(TN_pairs);
[h,stats_TS] = cdfplot(TS_pairs);
[h,stats_NN] = cdfplot(NN_pairs);
[h,stats_NS] = cdfplot(NS_pairs);
legend('Trained TrialNormal','Trained TrialShuffle','Novel TrialNormal','Novel TrialShuffle')



%% maintaining time series per bin as the multivariate... test across different bin sizes?


% need tools to deal with singularity of covariance matrix when n samples
% is smaller than p dimensions of the space (time samples or units for
% example)
%%%%% SOLUTION: shrinkage for sample covariance estimation

% for a subset of clsuters, 
% eliminate first 1/6 of response (1sec in good expts, 1.x in expts with
% mis-sampled rate
%for one subset of stimuli...

% calc mean and variance per stimulus
sigon = min(find(xtime_spikes>0));
sigoff = max(find(xtime_spikes<6));

calcstart = sigon + round((sigoff-sigon)/6); %size(TrialAvg_Stimulus,3)/6; %1/6 of way into 6-motif pseudo song

stimsize = floor(0.04*11025); %40msec

bins = [calcstart:stimsize:sigoff];

load('PopSumData_Normal_site1.mat')
ResponseSumFiltR_Normal = ResponseSumFiltR;

%novel songs
data = [];
stimind = 1;
for isong = 1:4
    thissong = squeeze(ResponseSumFiltR_Normal(isong,:,:));
    for ibin = 1:size(bins,2)-1
        data(stimind,:,:) = thissong(:,bins(ibin):(bins(ibin+1)-1));
        stimind = stimind+1;
    end
end

[m, v, s] = KLprep_Multivariate(data);
kl_pairwise_NN = KLpairwise_Multivariate(m , v);

%trained songs
data = [];
stimind = 1;
for isong = 5:8
    thissong = squeeze(ResponseSumFiltR_Normal(isong,:,:));
    for ibin = 1:size(bins,2)-1
        data(stimind,:,:) = thissong(:,bins(ibin):(bins(ibin+1)-1));
        stimind = stimind+1;
    end
end

[m, v] = KLprep_Multivariate(data);
kl_pairwise_TN = KLpairwise_Multivariate(m , v);
% figure;imagesc(KL_NetPseudoSyn_Trained_TrialNormal)
% title('KL_NetPseudoSyn_Trained_TrialNormal','Interpreter','none')


load('PopSumData_Shuffle_site1.mat')
ResponseSumFiltR_Shuffle = ResponseSumFiltR;

%novel songs
data = [];
stimind = 1;
for isong = 1:4
    thissong = squeeze(ResponseSumFiltR_Shuffle(isong,:,:));
    for ibin = 1:size(bins,2)-1
        data(stimind,:,:) = thissong(:,bins(ibin):(bins(ibin+1)-1));
        stimind = stimind+1;
    end
end

[m, v] = KLprep_Multivariate(data);
kl_pairwise_NS = KLpairwise_Multivariate(m , v);

%trained songs
data = [];
stimind = 1;
for isong = 5:8
    thissong = squeeze(ResponseSumFiltR_Shuffle(isong,:,:));
    for ibin = 1:size(bins,2)-1
        data(stimind,:,:) = thissong(:,bins(ibin):(bins(ibin+1)-1));
        stimind = stimind+1;
    end
end

[m, v] = KLprep_Univariate(data);
kl_pairwise_TS = KLpairwise_Univariate(m , v);
% figure;imagesc(KL_NetPseudoSyn_Trained_TrialNormal)
% title('KL_NetPseudoSyn_Trained_TrialNormal','Interpreter','none')


nsubsamps = 1000;
nreps = 1000;

NN_sumKL = KLsum(kl_pairwise_NN,nsubsamps,nreps);
TN_sumKL = KLsum(kl_pairwise_TN,nsubsamps,nreps);
NS_sumKL = KLsum(kl_pairwise_NS,nsubsamps,nreps);
TS_sumKL = KLsum(kl_pairwise_TS,nsubsamps,nreps);

figure;hold on
[h,stats_TN] = cdfplot(TN_sumKL);
[h,stats_TS] = cdfplot(TS_sumKL);
[h,stats_NN] = cdfplot(NN_sumKL);
[h,stats_NS] = cdfplot(NS_sumKL);
legend('Trained TrialNormal','Trained TrialShuffle','Novel TrialNormal','Novel TrialShuffle')
xlabel(['sum multivar KLdivergence for subset ' num2str(nsubsamps) ' stims; run ' num2st 'times on net pseudo-syn resp'])

%% power spectral density (welch)

load('metatoes.mat','xtime_spikes');
load('PopSumData_Normal_site1.mat');
ResponseSumFiltR_Normal = ResponseSumFiltR;

% calc mean and variance per stimulus
sigon = min(find(xtime_spikes>0));
sigoff = max(find(xtime_spikes<6));

calcstart = sigon + round((sigoff-sigon)/6); %size(TrialAvg_Stimulus,3)/6; %1/6 of way into 6-motif pseudo song

nstims = size(ResponseSumFiltR_Normal,1);
ntrials = size(ResponseSumFiltR_Normal,2);

for istim = 1:nstims
    
    thisstim = squeeze(ResponseSumFiltR_Normal(istim,:,calcstart:sigoff));
    thisstim = thisstim - mean(mean(thisstim,2));
    thisstim_avg = mean(thisstim,1);
    thisstim_resid = thisstim - repmat(thisstim_avg,ntrials,1);
    
    
    for itrial = 1:ntrials;
        
         figure;hold on
         
        [pxx,f] = pwelch(thistrial(itrial,:),[],[],[],11025);
        line(f,10*log10(pxx),'color','k')

        [pxx,f] = pwelch(thisstim_avg,[],[],[],11025);
        line(f,10*log10(pxx),'color','r')
        
        [pxx,f] = pwelch(thisstim_resid(itrial,:),[],[],[],11025);
        line(f,10*log10(pxx),'color','b')
        
        xlabel('Frequency (Hz)')
        ylabel('Magnitude (dB)')
    end
    
    
end

%%
load('KL_MultiVar_netPseudoSynResponse.mat')

%% Network Diameters

kl_NN = triu(ones(496), 1);
kl_NN(kl_NN==1) = kl_pairwise_NN;
kl_NN = (kl_NN + kl_NN');

kl_TN = triu(ones(496), 1);
kl_TN(kl_TN==1) = kl_pairwise_TN;
kl_TN = (kl_TN + kl_TN');

kl_TS = triu(ones(496), 1);
kl_TS(kl_TS==1) = kl_pairwise_TS;
kl_TS = (kl_TS + kl_TS');

kl_NS = triu(ones(496), 1);
kl_NS(kl_NS==1) = kl_pairwise_NS;
kl_NS = (kl_NS + kl_NS');

%%
adj = kl_TS;
diam=0;
for i=1:size(adj,1)
    
    n=length(adj);
    d = inf*ones(1,n); % distance s-all nodes
    d(i) = 0;    % s-s distance
    T = 1:n;    % node set with shortest paths not found

    while not(isempty(T))
     [dmin,ind] = min(d(T));
        for j=1:length(T)
            if adj(T(ind),T(j))>0 & d(T(j))>d(T(ind))+adj(T(ind),T(j))
                d(T(j))=d(T(ind))+adj(T(ind),T(j));
            end
        end 
        T = setdiff(T,T(ind));
    end
   
    diam = max([max(d),diam]);
end







