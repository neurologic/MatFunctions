cd('/mnt/cube/Ice/kperks/B1087/MetaMat/')
currentpath = pwd

load('IntracellularData.mat','IntracellularData')
dt = IntracellularData{1}.dt;
IntracellStimNames = IntracellularData{1}.stimnames;

load('metatoes.mat','allchanxz','metatoes','allclusterid')
% load('behavior.mat')
ExtracellStimNames = {};
for istim = 1:size(metatoes{1}.stims,1)
    ExtracellStimNames(istim) = metatoes{1}.stims{istim}.name;
end

load(['MVRegressPrep.mat'],'SignalFiltR','include_inds','isNarrow')
loc = allchanxz(include_inds);
allinds = [1:size(include_inds,1)];

CellID = 2; %which intracellular cell

Pop = SignalFiltR;
calcstart = size(Pop,3)/6;
nstims = size(SignalFiltR,2);
nclusters = size(SignalFiltR,1);

stiminds = [1:nstims];

subset_win = [calcstart+1:(2*calcstart)];
% subset_win = [calcstart+1:size(Pop,3)];

nsamps = size(subset_win,2);

Vm = [];
Meta_VarVm = [];

MSE = [];
Meta_Criteria_B = [];
Meta_B = [];
Meta_MeanMSE = [];
Meta_DF = [];
Meta_Vm = [];
VmData = [];
Meta_VmHat_Alpha = [];

%get global offset and norm for Vm across all stims
trialavgs = cellfun(@(x) mean(x,1),IntracellularData{CellID}.stim_data,'UniformOutput',0);
trialsamps =[];
for i = 1:size(trialavgs,2)
    trialsamps(i) = size(trialavgs{i},2) ;
end
for i = 1:size(trialavgs,2)
    IntracellularData{CellID}.stim_data{i} = IntracellularData{CellID}.stim_data{i}(:,1:min(trialsamps));
end
trialavgs = cellfun(@(x) mean(x,1),IntracellularData{CellID}.stim_data,'UniformOutput',0);

trialavgs = cell2mat(trialavgs');
globaloffset = min(min(trialavgs));
trialavgs = trialavgs - globaloffset;
globalnormfactor = max(max(trialavgs));
trialavgs = trialavgs ./ globalnormfactor;
Meta_Vm = trialavgs;
%Meta_Vm are the normalized traces that will go into regression calculation
%to maintain global tuning

MSE = [];
motifind = 1;
for istim = 1:nstims
    thisVm = IntracellularData{CellID}.stim_data{istim}(:,1:min(trialsamps));
    
    for imot = 1:5
        subset_win = [(imot*calcstart)+1:((imot+1)*calcstart)];
        
        thisMot = thisVm(:,subset_win);
        
        VmTrials = zscore(thisMot')';
%         
        VmAvg = [];
        nstrap = 50;
        ntrials = size(thisMot,1);
        nleavein = floor(0.75*ntrials);
        for istrap = 1:nstrap
            v = randperm(ntrials);
            VmAvg(istrap,:) = zscore(mean(thisMot(v(1:nleavein),:),1)')';            
        end
        
%         get pairwise mean squared error between bootstrap trial avg pairs
        strapind = 1;
        for istrap1 = 1:nstrap
            for istrap2 = istrap1+1:nstrap
                MSE{motifind}(strapind) = mean((VmAvg(istrap1,:) - VmAvg(istrap2,:)).^2);
                strapind = strapind+1;
                %             line([1:calcstart]*dt,((VmTrials(itrial,:) - VmAvg).^2)','color','k')
            end
        end
        
        figure; hold on
        line([1:size(VmAvg,2)]*dt,VmAvg(1,:)','color',[0.5,0.5,0.5])
        line([1:size(VmAvg,2)]*dt,mean(VmAvg,1)','color','r','LineWidth',2)
        legend('75% subset trial averages (n = 50)','trial average')
        line([1:size(VmAvg,2)]*dt,VmAvg','color',[0.5,0.5,0.5])
        line([1:size(VmAvg,2)]*dt,mean(VmAvg,1)','color','r','LineWidth',2)
        ylabel('zscore Vm')
        xlabel('seconds')
        title(['mean pairwise estimate MSE ' num2str(mean(MSE{motifind}))])


%         for istrap = 1:nstrap
%             MSE{motifind}(istrap) = mean((VmAvg(istrap,:) - zscore(mean(thisMot,1))).^2);
%         end
        
%         Vm MSE for goodness of CV for lasso
        
%         for itrial = 1:ntrials
%             MSE{motifind}(itrial) = mean((VmTrials(itrial,:) - zscore(mean(thisMot))).^2);
%             %             line([1:calcstart]*dt,((VmTrials(itrial,:) - VmAvg).^2)','color','k')
%         end
        %
        %         figure;hold on
        %         for itrial = 1:size(thisVm,1)
        %             MSE{motifind}(itrial) = mean((thisMot(itrial,:) - mean(thisMot)).^2);
        % %             line([1:calcstart]*dt,((thisMot(itrial,:) -  mean(thisMot)).^2)','color','k')
        %
        %         end
        %
        %
        
%         Meta_ZVm(:,motifind) = VmAvg';
        
        Meta_VarVm(motifind) = var(mean(thisMot));

        % sort the non-zscored but globally normalized Vm for regression
        VmData(motifind,:) = Meta_Vm(istim,subset_win);
        
        motifind = motifind+1;
    end
end
MeanMSE = cellfun(@(x) median(x),MSE)';

AllFits = [];
Meta_Criteria_B = [];
nmotifs = size(VmData,1);
ScaleMSE = 3;

for trialshuffleid = 1:10;
    % load('MVRegressPrep_Shuffled.mat','SignalFiltR','include_inds','isNarrow','ExtracellStimNames','IntracellStimNames')
    %     load('MVRegressPrep.mat','SignalFiltR','include_inds','isNarrow','ExtracellStimNames','IntracellStimNames')
    load(['MVRegressPrep_Subsample10Trial_' num2str(trialshuffleid) '.mat'],'SignalFiltR','include_inds','isNarrow')
    
    fileID = fopen('matlog.txt','w');
    t = datestr(datetime('now'));
    fprintf(fileID,'%s %i %s\n','on shuffle rep', trialshuffleid, t);
    fclose(fileID);
    Pop = SignalFiltR;
    
    %     Vm = zeros(nsamps,nstims);
    BetaWeightMatrix = zeros(nclusters,nmotifs);
    PVmData = zeros(nmotifs,nsamps,nclusters);
    VmHat = zeros(nmotifs,nsamps,nclusters);
    
    
    % get global normalization factor for each cluster across stims
    % (already offset to 0)
    % do regression off of globally normalized unit responses, not z-scored
    % to compare to globally normalized Vm
    tmpmax = abs((max(max(Pop,[],3),[],2)));
    tmpmin = abs((min(min(Pop,[],3),[],2)));
    cluster_globalnorm = max([tmpmax,tmpmin],[],2);
    
    PopNorm = [];
    for iunit = 1:size(Pop,1)
        PopNorm(iunit,:,:) = Pop(iunit,:,:) ./ cluster_globalnorm(iunit);
    end
    
    motifind = 1;
    PVmData = [];
   
    for istim = 1:nstims
        
        for imot = 1:5
            subset_win = [(imot*calcstart)+1:((imot+1)*calcstart)];
            % %
            Pop_stim = squeeze(PopNorm(:,stiminds(istim),subset_win))';
            Pop_stim(isnan(Pop_stim)) = 0;

            PVmData(motifind,:,:) = Pop_stim;
            
            motifind = motifind+1;
        end
    end
    %can sum across PVmData 2nd dimension to get some sense of lifetime sparseness of
    %each unit
%     tmp = squeeze(mean(PVmData,2));
%     tmp(find(tmp>0.1)) = 1;
%     figure;imagesc(tmp)
% %     colormap('gray')
%     title('binarized responses if mean>0.1')
%     ylabel('motif')
%     xlabel('unit')
%     
    VmHat_Alpha = [];
    mse_vm_ind = [];
    
    cv = 2;
    alphaVec = [1,0.01];
    %     lambdaVec = logspace(1,-4,50);
    for istim = 1:nmotifs
        
        thisVm = VmData(istim,:);

        thisPop = squeeze(PVmData(istim,:,:));
        %remove nearly zero psth that are causing empty Fit
%         thisPop(:,find(sum(thisPop,1)<0.00001)) = 0;
        thisPop(find(thisPop<0.0001)) = 0;
        
        for ialpha = [1,2];
            [B,FitInfo] = lasso(thisPop,thisVm,'Alpha',alphaVec(ialpha));% %now that don't use calculated MSE, maybe i can run this without CV and maybe it is faster then,'CV',cv);%,'Lambda',lambdaVec);
            nLambda = size(B,2);
            AllFits{trialshuffleid,ialpha,istim} = FitInfo;
            % lassoPlot(B,FitInfo,'PlotType','CV');
            % lassoPlot(B,FitInfo,'PlotType','Lambda');
            
            Meta_B(trialshuffleid,ialpha,istim,:,:) = B;

            
            %so that do not get an empty prediction
            %             if max(find(FitInfo.MSE <= MeanMSE(istim))) ==100
            %                 mse_vm_ind(ialpha,istim) = 99;
            %             end
            thisMSe = [];
            for ilambda = 1:nLambda
                VmHatdata = [];
                for iterm = 1:size(B,1)
                    VmHatdata(:,iterm) = B(iterm,ilambda)*thisPop(:,iterm);
                end
                thisMSE(ilambda) = mean((zscore(thisVm)' - zscore(sum(VmHatdata,2))).^2);
            end
            
%             mse_vm_ind(ialpha,istim) = max(find(thisMSE <= MeanMSE(istim)));
%             mse_vm_ind(ialpha,istim) = max(find(thisMSE <= (min(thisMSE)+(2*MeanMSE(istim)))));
            
            mse_vm_ind(ialpha,istim) = max(find(thisMSE <= (min(thisMSE)+(ScaleMSE*MeanMSE(istim)))));
            
            Meta_DF(trialshuffleid,ialpha,istim) = FitInfo.DF(mse_vm_ind(ialpha,istim));
            
            Meta_Criteria_B(trialshuffleid,ialpha,istim,:) = B(:,mse_vm_ind(ialpha,istim));
            VmHatdata = [];
            for iterm = 1:size(B,1)
                VmHatdata(:,iterm) = B(iterm,mse_vm_ind(ialpha,istim))*thisPop(:,iterm);
            end
            
            VmHat_Alpha(ialpha,:,istim) = sum(VmHatdata,2);
            
            
        end
    end
    
    Meta_Criteria_B = Meta_Criteria_B;
    Meta_B = Meta_B;
    Meta_DF = Meta_DF;
    Meta_PVm(trialshuffleid,:,:,:) = PVmData;
    Meta_VmHat_Alpha(trialshuffleid,:,:,:) = VmHat_Alpha;
%     Meta_VmHat_Ridge(trialshuffleid,:,:) = squeeze(sum(VmHat_Ridge,3));
    
end

Meta_Vm = Vm;
Meta_MeanMSE = MeanMSE;
Meta_VarVm = Meta_VarVm;

% save(['Lasso_ClusterID_Overlap_MetaVectors_CellID' num2str(CellID) '.mat'],'Meta_B','Meta_Criteria_B','Meta_MeanMSE','MetaPVm','Vm','AllFits','alphaVec','lambdaVec')
save(['Lasso_ClusterID_Overlap_MetaVectors_CellID' num2str(CellID) '_3MSECriteria.mat'],'dt','alphaVec','Meta_VarVm','Meta_B','Meta_Criteria_B','Meta_MeanMSE','Meta_Vm','VmData','Meta_VmHat_Alpha','Meta_DF','AllFits')