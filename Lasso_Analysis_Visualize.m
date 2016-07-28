
BirdID = {'B1075'
    'B1083'
    'B1083'
    'B1087'
    'B1087'
    'B1112'
    'B1112'
    'B1112'};

Sample = {'Lasso_ClusterID_Overlap_MetaVectors_CellID1_3MSECriteria.mat'
    'Lasso_ClusterID_Overlap_MetaVectors_CellID1_3MSECriteria.mat'
    'Lasso_ClusterID_Overlap_MetaVectors_CellID2_3MSECriteria.mat'
    'Lasso_ClusterID_Overlap_MetaVectors_CellID1_3MSECriteria.mat'
    'Lasso_ClusterID_Overlap_MetaVectors_CellID2_3MSECriteria.mat'
    'Lasso_ClusterID_Overlap_MetaVectors_CellID1_3MSECriteria.mat'
    'Lasso_ClusterID_Overlap_MetaVectors_CellID2_3MSECriteria.mat'
    'Lasso_ClusterID_Overlap_MetaVectors_CellID3_3MSECriteria.mat'};
%%
Lalpha = 1;
Ralpha = 2;

DF_Ratio = [];
MCorr_lasso_stim_ALL = [];
MCorr_lasso_subsample_ALL = [];
MCorr_ridge_stim_ALL = [];
MCorr_ridge_subsample_ALL = [];
    
for pairID = 1:size(BirdID,1);
    cd(['/Users/kperks/mnt/cube/Ice/kperks/' BirdID{pairID} '/MetaMat/'])
    load(Sample{pairID})
    
    nsubsample = size(Meta_VmHat_Alpha,1);
    nstims = size(Meta_VmHat_Alpha,4);
    
    % Lasso DF to Ridge DF ( / total DF?)
    
    DF_Ratio{pairID} = squeeze(mean(Meta_DF(:,Lalpha,:),1)) ./ squeeze(mean(Meta_DF(:,Ralpha,:),1));
    
    % SIMILARITY
    
    edges = [-1:0.05:1];
    
    %%%%LASSO
    %Across stimuli
    MCorr_lasso_stim = [];
    for ishuffle = 1:nsubsample
        B = squeeze(Meta_Criteria_B(ishuffle,Lalpha,:,:));
        B_binary = B;
        B_binary(find(B_binary)) = 1;
        
        pairind = 1;
        for istim1 = 1:nstims -1
            for istim2 = istim1 +1 :nstims
                R1 = matthewscorr(B_binary(istim1,:), B_binary(istim2,:));
                R2 = matthewscorr(B_binary(istim2,:), B_binary(istim1,:));
                MCorr_lasso_stim(ishuffle,pairind) = mean([R1,R2]);
                
                pairind = pairind+1;
            end
        end
    end
    % n_MCC_lasso_stim = histc(mean(MCorr_stim,1),edges) ./ size(MCorr_stim,2);
    
    %Across subsampled datasets
    MCorr_lasso_subsample = [];
    for istim = 1:nstims
        B = squeeze(Meta_Criteria_B(:,Lalpha,istim,:));
        B_binary = B;
        B_binary(find(B_binary)) = 1;
        
        shufflepair= 1;
        for ishuffle1 = 1:nsubsample-1
            for ishuffle2 = ishuffle1 +1 :nsubsample
                R1 = matthewscorr(B_binary(ishuffle1,:), B_binary(ishuffle2,:));
                R2 = matthewscorr(B_binary(ishuffle2,:), B_binary(ishuffle1,:));
                MCorr_lasso_subsample(istim,shufflepair) = mean([R1,R2]);
                
                shufflepair = shufflepair+1;
            end
        end
    end
    % n_MCC_lasso_subsample = histc(mean(MCorr_subsample,1),edges) ./ size(MCorr_subsample,2);
    
    normfactor = max(mean(MCorr_lasso_subsample,1));
    MCorr_lasso_stim_ALL{pairID} = mean(MCorr_lasso_stim,1) ./ normfactor;
    MCorr_lasso_subsample_ALL{pairID} = mean(MCorr_lasso_subsample,1) ./ normfactor;
    
    %%%%RIDGE
    %Across stimuli
    MCorr_ridge_stim = [];
    for ishuffle = 1:nsubsample
        B = squeeze(Meta_Criteria_B(ishuffle,Ralpha,:,:));
        B_binary = B;
        B_binary(find(B_binary)) = 1;
        
        pairind = 1;
        for istim1 = 1:nstims -1
            for istim2 = istim1 +1 :nstims
                R1 = matthewscorr(B_binary(istim1,:), B_binary(istim2,:));
                R2 = matthewscorr(B_binary(istim2,:), B_binary(istim1,:));
                MCorr_ridge_stim(ishuffle,pairind) = mean([R1,R2]);
                
                pairind = pairind+1;
            end
        end
    end
    % n_MCC_ridge_stim = histc(mean(MCorr_stim,1),edges) ./ size(MCorr_stim,2);
    
    %Across subsampled datasets
    MCorr_ridge_subsample = [];
    for istim = 1:nstims
        B = squeeze(Meta_Criteria_B(:,Ralpha,istim,:));
        B_binary = B;
        B_binary(find(B_binary)) = 1;
        
        shufflepair= 1;
        for ishuffle1 = 1:nsubsample-1
            for ishuffle2 = ishuffle1 +1 :nsubsample
                R1 = matthewscorr(B_binary(ishuffle1,:), B_binary(ishuffle2,:));
                R2 = matthewscorr(B_binary(ishuffle2,:), B_binary(ishuffle1,:));
                MCorr_ridge_subsample(istim,shufflepair) = mean([R1,R2]);
                
                shufflepair = shufflepair+1;
            end
        end
    end
    % n_MCC_ridge_subsample = histc(mean(MCorr_subsample,1),edges) ./ size(MCorr_subsample,2);
    
    normfactor = max(mean(MCorr_ridge_subsample,1));
    MCorr_ridge_stim_ALL{pairID} = mean(MCorr_ridge_stim,1) ./ normfactor;
    MCorr_ridge_subsample_ALL{pairID} = mean(MCorr_ridge_subsample,1) ./ normfactor;
    

end


save('MetaLasso_AllPairs.mat', 'DF_Ratio', 'MCorr_lasso_stim_ALL', 'MCorr_lasso_subsample_ALL', 'MCorr_ridge_stim_ALL', 'MCorr_ridge_subsample_ALL');

%DF SUMMS AND PLOTS
edges = [0:0.05:1];

for ipair = 1:size(DF_Ratio,2)
    n(ipair,:) = histc(DF_Ratio{ipair},edges) ./ size(DF_Ratio{ipair},1);

end
figure;hold on
stairs(edges, mean(n,1),'LineWidth',4)
median(DF_Ratio


%SIMILARITY SUMMS AND PLOTS
median(MCorr_ridge_stim)

n_normMCC_lasso_stim = histc(MCorr_lasso_stim_ALL,edges) ./ size(MCorr_lasso_stim_ALL,2);
n_normMCC_lasso_subsample = histc(MCorr_lasso_subsample_ALL,edges) ./ size(MCorr_lasso_subsample_ALL,2);
n_normMCC_ridge_stim = histc(MCorr_ridge_stim_ALL,edges) ./ size(MCorr_ridge_stim_ALL,2);
n_normMCC_ridge_subsample = histc(MCorr_ridge_subsample_ALL,edges) ./ size(MCorr_ridge_subsample_ALL,2);

figure;hold on
stairs(edges,n_normMCC_ridge_subsample,'color','k','LineWidth',3)
stairs(edges,n_normMCC_ridge_stim,'color','b','LineWidth',3)
title('Ridge Regularization')
xlabel('MCC normalized to max within dataset')
ylabel('p(MCC)')
legend('Stimulus Condition','Subsampled Condition')
axis tight
ylims = get(gca,'YLim');
for n = 1:size(MCorr_subsample,2)
    scatter(median(mean(MCorr_subsample,1)),max(ylims)+0.02,200,'*','k')
    scatter(median(mean(MCorr_stim,1)),max(ylims)+0.02,200,'*','b')
end
set(gca,'YLim',[0,max(ylims)+0.1])

figure;hold on
stairs(edges,n_normMCC_lasso_subsample,'color','k','LineWidth',3)
stairs(edges,n_normMCC_lasso_stim,'color','b','LineWidth',3)
title('Lasso Regularization')
xlabel('MCC normalized to max within dataset')
ylabel('p(MCC)')
legend('Stimulus Condition','Subsampled Condition')
axis tight
ylims = get(gca,'YLim');
for n = 1:size(MCorr_lasso_subsample,2)
    scatter(median(mean(MCorr_lasso_subsample,1)),max(ylims)+0.02,200,'*','k')
    scatter(median(mean(MCorr_lasso_stim,1)),max(ylims)+0.02,200,'*','b')
end
set(gca,'YLim',[0,max(ylims)+0.1])


%% FIGURES for Lasso DF to Ridge DF ( / total DF?)
figure;
hold on
scatter(Meta_MeanMSE,squeeze(mean(Meta_DF(:,Lalpha,:),1)),200,'b')
scatter(Meta_MeanMSE,squeeze(mean(Meta_DF(:,Ralpha,:),1)),200,'r')
ylabel('degrees of freedom')
xlabel('MSE(Vm pairwise trials)')

figure;
hold on
scatter(Meta_VarVm,squeeze(mean(Meta_DF(:,Lalpha,:),1)),200,'b')
scatter(Meta_VarVm,squeeze(mean(Meta_DF(:,Ralpha,:),1)),200,'r')
ylabel('degrees of freedom')
xlabel('Var(Vm temporal)')

figure;
hold on
scatter(Meta_MeanMSE./Meta_VarVm',squeeze(mean(Meta_DF(:,Lalpha,:),1)),200,'b')
scatter(Meta_MeanMSE./Meta_VarVm',squeeze(mean(Meta_DF(:,Ralpha,:),1)),200,'r')
ylabel('degrees of freedom')
xlabel('MSE(Vm pairwise trials) / Var(Vm temporal)')

figure;
hold on
scatter(Meta_MeanMSE./Meta_VarVm',squeeze(mean(Meta_DF(:,Lalpha,:),1))./squeeze(mean(Meta_DF(:,Ralpha,:),1)),200,'k','fill')
ylabel('DF under lasso / DF under ridge')
xlabel('MSE (Vm pairwise trials) / Var (Vm temporal)')

figure;
scatter(squeeze(mean(Meta_DF(:,Lalpha,:),1)),squeeze(mean(Meta_DF(:,Ralpha,:),1)),200,'k','fill')
xlabel('DF under Lasso')
ylabel('DF under Ridge')
axis tight
SetAxisUnity(gcf)

figure;scatter(Meta_MeanMSE,Meta_VarVm,200,'k','fill')
xlabel('MSE')
ylabel('Variance Temporal')

%%
imot  =1;
%
% ialpha = 2;
% figure;plot(squeeze(mean(Meta_VmHat_Alpha(:,ialpha,:,imot),1)))
% title('mean prediction (across shuffles) for Motif 1 under Ridge')
%
% ialpha = 1;
% figure;plot(squeeze(mean(Meta_VmHat_Alpha(:,ialpha,:,imot),1)))
% title('mean prediction (across shuffles) for Motif 1 under Lasso')
%
% figure;
% plot(squeeze(VmData(imot,:)))
% title('mean Vm for Motif 1 intracellular')

xtime = [1:size(VmData,2)].*dt;
figure; hold on
v = squeeze(VmData(imot,:));
v = zscore(v);
line(xtime,v,'color','k')
ialpha = 1;
l = squeeze(mean(Meta_VmHat_Alpha(:,ialpha,:,imot),1));
l = zscore(l);
line(xtime,l,'color','b')
ialpha = 2;
r = squeeze(mean(Meta_VmHat_Alpha(:,ialpha,:,imot),1));
r = zscore(r);
line(xtime,r,'color','r')
legend('Vm',['Lasso DF ' num2str(mean(Meta_DF(:,1,imot)))],['Ridge DF ' num2str(mean(Meta_DF(:,2,imot)))])
ylabel('zscore')
xlabel('sec')
title(['motif#' num2str(imot)])


%% no longer doing CV version of the analysis so don't get these plots.... can still get the lambda
ialpha = 1;
trialshuffleid = 1;
imot = 1
FitInfo = AllFits{trialshuffleid,ialpha,imot} ;
B = squeeze(Meta_B(trialshuffleid,ialpha,imot,:,:));
lassoPlot(B,FitInfo,'PlotType','Lambda');

ialpha = 2;
trialshuffleid = 1;
imot = 1
FitInfo = AllFits{trialshuffleid,ialpha,imot} ;
B = squeeze(Meta_B(trialshuffleid,ialpha,imot,:,:));
lassoPlot(B,FitInfo,'PlotType','Lambda');

%%
% for each shuffle rep get Beta vector at criteria
% average across shuffle conditions to get average Beta Vector
% ialpha = 1;
% B = squeeze(mean(Meta_Criteria_B(:,ialpha,:,:),1));
% B_binary = B;
% B_binary(find(B_binary)) = 1;
% figure;
% imagesc(B_binary)
% title('Lasso')
% ylabel('motif ID')
% xlabel('unit ID')
%
% ialpha = 2;
% B = squeeze(mean(Meta_Criteria_B(:,ialpha,:,:),1));
% B_binary = B;
% B_binary(find(B_binary)) = 1;
% figure;
% imagesc(B_binary)
% title('Ridge')
% ylabel('motif ID')
% xlabel('unit ID')

%%%%% to plot can I sort these by location? (might want to organize it
%%%%% this way anyway for next analysis)
load(['MVRegressPrep.mat'],'include_inds')
load('metatoes.mat','allchanxz')
loc = allchanxz(include_inds);

[sorted_location,isort] = sort(loc);

ialpha = 1;
B = squeeze(mean(Meta_Criteria_B(1,ialpha,:,:),1));
B_binary = B;
B_binary(find(B_binary)) = 1;
figure;
B_sorted = B_binary(:,isort);
imagesc(B_sorted)
ticks = get(gca,'XTick');
title('Lasso')
ylabel('motif ID')
xlabel('units by location ID')
set(gca,'XTickLabel',sorted_location(ticks))

ialpha = 2;
B = squeeze(mean(Meta_Criteria_B(1,ialpha,:,:),1));
B_binary = B;
B_binary(find(B_binary)) = 1;
figure;
B_sorted = B_binary(:,isort);
imagesc(B_sorted)
ticks = get(gca,'XTick');
title('Ridge')
ylabel('motif ID')
xlabel('units by location ID')
set(gca,'XTickLabel',sorted_location(ticks))



%%%% Analogue version is not clearly visualized for the mean across
%%%% subsampled datasets
% ialpha = 1;
% B = squeeze(mean(Meta_Criteria_B(:,ialpha,:,:),1));
% figure;
% B_sorted = B(:,isort);
% imagesc(B_sorted)
% ticks = get(gca,'XTick');
% title('Lasso')
% ylabel('motif ID')
% xlabel('units by location ID')
% set(gca,'XTickLabel',sorted_location(ticks))
%
% ialpha = 2;
% B = squeeze(mean(Meta_Criteria_B(1,ialpha,:,:),1));
% figure;
% B_sorted = B(:,isort);
% imagesc(B_sorted)
% ticks = get(gca,'XTick');
% title('Ridge')
% ylabel('motif ID')
% xlabel('units by location ID')
% set(gca,'XTickLabel',sorted_location(ticks))


%%
% compare similarity across motifs
%%%%%%%% similarity of mean Beta at Criteria
%%%%%%%% mean similarity of Beta at Criteria (similarity across motifs calc
%%%%%%%% under each shuffle and then averaged across shuffles at the end

%%%%%%%%%%%%%%%%%%%%
%%%%% SIMILARITY BETWEEN STIMULI (averaged across shuffles)
%%%%%%%%%%%%%%%%%%%%
edges = [-1:0.05:1];

idx = true(2);
idx = ~tril(idx);

ialpha = 1; %for lasso
MCorr = [];
Corrcoef = [];
for ishuffle = 1:nsubsample
    B = squeeze(Meta_Criteria_B(ishuffle,ialpha,:,:));
    B_binary = B;
    B_binary(find(B_binary)) = 1;
    
    pairind = 1;
    for istim1 = 1:nstims -1
        for istim2 = istim1 +1 :nstims
            R1 = matthewscorr(B_binary(istim1,:), B_binary(istim2,:));
            R2 = matthewscorr(B_binary(istim2,:), B_binary(istim1,:));
            MCorr(ishuffle,pairind) = mean([R1,R2]);
            
            c = corrcoef(B(istim1,:), B(istim2,:));
            Corrcoef(ishuffle,pairind) = c(idx);
            
            pairind = pairind+1;
        end
    end
end
n_MCC_lasso = histc(mean(MCorr,1),edges) ./ size(MCorr,2);
n_CC_lasso = histc(mean(Corrcoef,1),edges) ./ size(Corrcoef,2);
figure;
hold on
stairs(edges,n_CC_lasso,'color','r')
stairs(edges,n_MCC_lasso,'color','b')
legend('Correlation Coefficients','Mathews')
title('lasso reg - mean pairwise motif similarity')

ialpha = 2; %for ridge
MCorr = [];
Corrcoef = [];
for ishuffle = 1:nsubsample
    B = squeeze(Meta_Criteria_B(ishuffle,ialpha,:,:));
    B_binary = B;
    B_binary(find(B_binary)) = 1;
    
    pairind = 1;
    for istim1 = 1:nstims -1
        for istim2 = istim1 +1 :nstims
            R1 = matthewscorr(B_binary(istim1,:), B_binary(istim2,:));
            R2 = matthewscorr(B_binary(istim2,:), B_binary(istim1,:));
            MCorr(ishuffle,pairind) = mean([R1,R2]);
            
            c = corrcoef(B(istim1,:), B(istim2,:));
            Corrcoef(ishuffle,pairind) = c(idx);
            
            pairind = pairind+1;
        end
    end
end
n_MCC_ridge = histc(mean(MCorr,1),edges) ./ size(MCorr,2);
n_CC_ridge = histc(mean(Corrcoef,1),edges) ./ size(Corrcoef,2);
figure;
hold on
stairs(edges,n_CC_ridge,'color','r')
stairs(edges,n_MCC_ridge,'color','b')
legend('Correlation Coefficients','Mathews')
title('ridge reg - mean pairwise motif similarity per subsample')

% figure;hold on; stairs(edges,n_MCC_lasso,'color','b'); stairs(edges,n_MCC_ridge,'color','r')
% figure;hold on; stairs(edges,n_CC_lasso,'color','b'); stairs(edges,n_CC_ridge,'color','r')


% compare similarity across shuffles for same motif
%%%%%%%%

ialpha = 1;
MCorr = [];
Corrcoef = [];
for istim = 1:nstims
    B = squeeze(Meta_Criteria_B(:,ialpha,istim,:));
    B_binary = B;
    B_binary(find(B_binary)) = 1;
    
    shufflepair= 1;
    for ishuffle1 = 1:nsubsample-1
        for ishuffle2 = ishuffle1 +1 :nsubsample
            R1 = matthewscorr(B_binary(ishuffle1,:), B_binary(ishuffle2,:));
            R2 = matthewscorr(B_binary(ishuffle2,:), B_binary(ishuffle1,:));
            MCorr(istim,shufflepair) = mean([R1,R2]);
            
            c = corrcoef(B(ishuffle1,:), B(ishuffle2,:));
            Corrcoef(istim,shufflepair) = c(idx);
            
            shufflepair = shufflepair+1;
        end
    end
end
n_MCC_lasso_shuffle = histc(mean(MCorr,1),edges) ./ size(MCorr,2);
n_CC_lasso_shuffle = histc(mean(Corrcoef,1),edges) ./ size(Corrcoef,2);
figure;hold on;
stairs(edges,n_MCC_lasso_shuffle,'color','c');
stairs(edges,n_CC_lasso_shuffle,'color','m')
legend('matthews','corrcoef')
title('lasso reg - mean pairwise subsample similarity per motif')

ialpha = 2;
MCorr = [];
Corrcoef = [];
for istim = 1:nstims
    B = squeeze(Meta_Criteria_B(:,ialpha,istim,:));
    B_binary = B;
    B_binary(find(B_binary)) = 1;
    
    shufflepair= 1;
    for ishuffle1 = 1:nsubsample-1
        for ishuffle2 = ishuffle1 +1 :nsubsample
            R1 = matthewscorr(B_binary(ishuffle1,:), B_binary(ishuffle2,:));
            R2 = matthewscorr(B_binary(ishuffle2,:), B_binary(ishuffle1,:));
            MCorr(istim,shufflepair) = mean([R1,R2]);
            
            c = corrcoef(B(ishuffle1,:), B(ishuffle2,:));
            Corrcoef(istim,shufflepair) = c(idx);
            
            shufflepair = shufflepair+1;
        end
    end
end
n_MCC_ridge_shuffle = histc(mean(MCorr,1),edges) ./ size(MCorr,2);
n_CC_ridge_shuffle = histc(mean(Corrcoef,1),edges) ./ size(Corrcoef,2);
figure;hold on;
stairs(edges,n_MCC_ridge_shuffle,'color','c');
stairs(edges,n_CC_ridge_shuffle,'color','m')
legend('matthews','corrcoef')
title('ridge reg - mean pairwise subsample similarity per motif')

figure;hold on;
stairs(edges,n_MCC_lasso,'color','b')
stairs(edges,n_MCC_lasso_shuffle,'color','k')
xlabel('MCC')
ylabel('p(MCC)')
title('lasso regularization')
legend('MCC across motifs','MCC across subsample bootstraps')

figure;hold on;
stairs(edges,n_MCC_ridge,'color','r')
stairs(edges,n_MCC_ridge_shuffle,'color','k')
xlabel('MCC')
ylabel('p(MCC)')
title('ridge regularization')
legend('MCC across motifs','MCC across subsample bootstraps')

figure;hold on;
stairs(edges,n_CC_lasso,'color','b')
stairs(edges,n_CC_lasso_shuffle,'color','k')
xlabel('CC')
ylabel('p(CC)')
title('lasso regularization')
legend('CC across motifs','CC across subsample bootstraps')

figure;hold on;
stairs(edges,n_CC_ridge,'color','r')
stairs(edges,n_CC_ridge_shuffle,'color','k')
xlabel('CC')
ylabel('p(CC)')
title('ridge regularization')
legend('CC across motifs','CC across subsample bootstraps')

%
% ialpha = 2;
% figure;
% imagesc(squeeze(Meta_Criteria_B(:,ialpha,istim,:)));
% ylabel('subset iteration')
% xlabel('unit ID')
% title('motif#1 analog B')

ialpha = 1;
figure;
B = squeeze(Meta_Criteria_B(:,ialpha,istim,:));
B_binary = B;
B_binary(find(B_binary)) = 1;
imagesc(B_binary);
ylabel('subset iteration')
xlabel('unit ID')
title('motif#1 binary B')


%%
%%%%%% taking mean vector B across subsample iterations before doing
%%%%%% similarity calculations...
idx = true(2);
idx = ~tril(idx);
edges = [-1:0.05:1];

ialpha = 1;
MCorr = [];
Corrcoef = [];
B = squeeze(mean(Meta_Criteria_B(:,ialpha,:,:),1));
B_binary = B;
B_binary(find(B_binary)) = 1;
pairind = 1;
for istim1 = 1:nstims -1
    for istim2 = istim1 +1 :nstims
        R1 = matthewscorr(B_binary(istim1,:), B_binary(istim2,:));
        R2 = matthewscorr(B_binary(istim2,:), B_binary(istim1,:));
        MCorr(pairind) = mean([R1,R2]);
        
        c = corrcoef(B(istim1,:), B(istim2,:));
        Corrcoef(pairind) = c(idx);
        
        pairind = pairind+1;
    end
end
n_MCC_lasso = histc(MCorr,edges) ./ size(MCorr,2);
n_CC_lasso = histc(Corrcoef,edges) ./ size(Corrcoef,2);
figure;
hold on
stairs(edges,n_CC_lasso,'color','r')
stairs(edges,n_MCC_lasso,'color','b')
legend('Correlation Coefficients','Mathews')
title('lasso regularization')

ialpha = 2;
MCorr = [];
Corrcoef = [];
B = squeeze(mean(Meta_Criteria_B(:,ialpha,:,:),1));
B_binary = B;
B_binary(find(B_binary)) = 1;
pairind = 1;
for istim1 = 1:nstims -1
    for istim2 = istim1 +1 :nstims
        R1 = matthewscorr(B_binary(istim1,:), B_binary(istim2,:));
        R2 = matthewscorr(B_binary(istim2,:), B_binary(istim1,:));
        MCorr(pairind) = mean([R1,R2]);
        
        c = corrcoef(B(istim1,:), B(istim2,:));
        Corrcoef(pairind) = c(idx);
        
        pairind = pairind+1;
    end
end
n_MCC_ridge = histc(MCorr,edges) ./ size(MCorr,2);
n_CC_ridge = histc(Corrcoef,edges) ./ size(Corrcoef,2);
figure;
hold on
stairs(edges,n_CC_ridge,'color','r')
stairs(edges,n_MCC_ridge,'color','b')
legend('Correlation Coefficients','Mathews')
title('ridge regularization')


figure;hold on; stairs(edges,n_MCC_lasso,'color','b'); stairs(edges,n_MCC_ridge,'color','r')
legend('lasso','ridge')
title('Matthews')
figure;hold on; stairs(edges,n_CC_lasso,'color','b'); stairs(edges,n_CC_ridge,'color','r')
legend('lasso','ridge')
title('Corrcoef')



%%%%%%%%%%%%%%%%%%%%
%%%%% SIMILARITY BETWEEN STIMULI (averaged across shuffles)
%%%%% BUT FIRST TAKING AVERAGE BETA VECTOR ACROSS SHUFFLES PER STIM
%%%%% BEFORE CORRELATE
%%%%%%%%%%%%%%%%%%%%
idx = true(2);
idx = ~tril(idx);

ialpha = 1;
LassoBCriteria_est = squeeze(sum(Meta_Criteria_B(:,ialpha,:,:),1));
ialpha = 2;
RidgeBCriteria_est = squeeze(sum(Meta_Criteria_B(:,ialpha,:,:),1));
pairind = 1;
for istim1 = 1:nstims -1
    for istim2 = istim1 +1 :nstims
        
        c = corrcoef(LassoBCriteria_est(istim1,:), LassoBCriteria_est(istim2,:));
        Corrcoef_Lasso_stim_est(pairind) = c(idx);
        
        c = corrcoef(RidgeBCriteria_est(istim1,:), RidgeBCriteria_est(istim2,:));
        Corrcoef_Ridge_stim_est(pairind) = c(idx);
        
        pairind = pairind+1;
    end
end






% for each motif across shuffles, get pairwise similarity... this is
% control distribution of dissimilarity expected given same motif, but different experiemnt


%% get spectral* entropy of spatial distribution beta weights
BirdID = {'B1075'
    'B1083'
    'B1083'
    'B1087'
    'B1087'
    'B1112'
    'B1112'
    'B1112'};

Sample = {'ClusterID_Overlap_MetaVectors_CellID1.mat'
    'Lasso_ClusterID_Overlap_MetaVectors_CellID1.mat'
    'Lasso_ClusterID_Overlap_MetaVectors_CellID2.mat'
    'Lasso_ClusterID_Overlap_MetaVectors_CellID1.mat'
    'Lasso_ClusterID_Overlap_MetaVectors_CellID2.mat'
    'Lasso_ClusterID_Overlap_MetaVectors_CellID1.mat'
    'Lasso_ClusterID_Overlap_MetaVectors_CellID2.mat'
    'Lasso_ClusterID_Overlap_MetaVectors_CellID3.mat'};
%%
for pairID = 2;%1:size(BirdID,1);
    cd(['/Users/kperks/mnt/cube/Ice/kperks/' BirdID{pairID} '/MetaMat/'])
    load(Sample{pairID})
    
    nsubsample = size(Meta_VmHat_Alpha,1);
    nstims = size(Meta_VmHat_Alpha,4);
    
    load('MVRegressPrep.mat','include_inds')
    load('metatoes.mat','allchanxz','allclusterid')
    
    these_locations = allchanxz(include_inds);
    locationVec = unique(these_locations);
    
    BetaVec = zeros(nsubsample,nstims,size(locationVec,1),1);
    BetaVec_shuffle = zeros(nsubsample,nstims,size(locationVec,1),1);
    
    Entropy = nan(nsubsample,nstims);
    Entropy_shuffle = nan(nsubsample,nstims);
    Entropy_ranked = nan(nsubsample,nstims);
    
    ialpha = 2; %for ridge condition
    for ishuffle = 1:nsubsample
        
        for istim = 1:nstims
            %         stimB = abs(squeeze(MetaBetaWeightMatrix(i,:,istim)));
            stimB = abs(squeeze(Meta_Criteria_B(ishuffle,ialpha,istim,:)));
            
            
            shuffleinds = randperm(size(stimB,1));
            stimB_shuffle = stimB(shuffleinds);
            %         [x,i] = sort(stimB);
            nloc = [];
            for l = 1:size(locationVec,1)
                linds = find(these_locations == locationVec(l));
                nloc(l) = size(linds,1);
                if ~isempty(linds)
                    BetaVec(ishuffle,istim,l) = mean(stimB(linds));
                    BetaVec_shuffle(ishuffle,istim,l) = mean(stimB_shuffle(linds));
                end
            end
            
            %%%%%Shuffle control
            d = squeeze(BetaVec_shuffle(ishuffle,istim,:));
            d=d/sum(d); %Normalization
            %Entropy Calculation
            logd = log2(d);
            logd(isinf(logd)) = 0;
            Entropy_shuffle(ishuffle,istim) = -sum(d.*logd)/log2(size(d,1));
            
            %%%%%Non-shuffle condition
            d = squeeze(BetaVec(ishuffle,istim,:));
            d=d/sum(d); %Normalization
            %Entropy Calculation
            logd = log2(d);
            logd(isinf(logd)) = 0;
            Entropy(ishuffle,istim) = -sum(d.*logd)/log2(size(d,1));
            
            %%%%%Rank Order Beta by location to get max peakiness
            sortlocN = nloc; %sort(nloc,'descend'); %
            sortB = sort(stimB,'descend');
            distB = [];
            for l = 1:size(locationVec,1)
                if ~isempty(find(these_locations == locationVec(l)))
                    distB(l) = mean(sortB(1:sortlocN(l)));
                    sortB = sortB(sortlocN(l)+1:end);
                end
            end
            d_rank = distB ./ sum(distB);%Normalization
            %Entropy Calculation
            logd = log2(d_rank);
            logd(isinf(logd)) = 0;
            Entropy_ranked(ishuffle,istim) = -sum(d_rank.*logd)/log2(size(d_rank,2));
            
        end
    end
    
    Entropy_meta{pairID} = mean(Entropy,1)';
    Entropy_meta_Shuffle{pairID} = mean(Entropy_shuffle,1)';
    Entropy_meta_Ranked{pairID} = mean(Entropy_ranked,1)';
end

figure;plot(locationVec,d)
figure;plot(locationVec,d_rank)
% after getting estimate entropy for all samples... and copy/paste from excel doc
%for each sample, plot relative to minimum possible entropy(most peaky possible)
e = Entropy_meta;
s = Entropy_meta_Shuffle;
r = Entropy_meta_Ranked;

alls = [];
alle = [];
allr = [];
for n = 1:size(s,2)
    alls = [alls;s{n}];
    alle = [alle;e{n}];
    allr = [allr;r{n}];
end
meanEntropy = mean(alle)
getCDFconf(alle,95)
meanEntropyShuffle = mean(alls)
getCDFconf(alls,95)
meanEntropyRanked = mean(allr)
getCDFconf(allr,95)
meanEntropy_RatioToRanked = mean(alle-allr)
getCDFconf(alle-allr,95)
meanEntropyShuffle_RatioToRanked = mean(alls-allr)
getCDFconf(alls-allr,95)

edges = [0:0.05:1];
nE = [];
nS = [];
nR = [];
for n = 1:size(e,2)
    nE(:,n) = histc(e{n},edges) ./size(e{n},1);
    nS(:,n) = histc(s{n},edges) ./size(s{n},1);
    nR(:,n) = histc(r{n},edges) ./size(r{n},1);
end
figure;hold on
stairs(edges,mean(nR,2),'color','r','LineWidth',4)
stairs(edges,mean(nS,2),'color','b','LineWidth',4)
stairs(edges,mean(nE,2),'color','k','LineWidth',4)
xlabel('Normalized Entropy')
ylabel('probability(Entropy)')
legend('Maximally Peaked','Shuffle Condition','Non-Shuffle Condition')
axis tight
ylims = get(gca,'YLim');
% set(gca,'YLim',[0,0.15])
set(gca,'YLim',[0,ylims(2)+0.05])
for n = 1:size(s,2)
    scatter(mean(e{n}),max(ylims)+0.02,200,'*','b')
    scatter(mean(s{n}),max(ylims)+0.02,200,'*','k')
    scatter(mean(r{n}),max(ylims)+0.02,200,'*','r')
end

edges = [0:0.05:1];
nEvR = [];
nSvR = [];
for n = 1:size(r,2)
    nEvR(:,n) = histc(r{n}./e{n},edges) ./size(r{n},1);
    nSvR(:,n) = histc(r{n}./s{n},edges) ./size(r{n},1);
end
figure;hold on
stairs(edges,mean(nEvR,2),'color','b','LineWidth',4)
stairs(edges,mean(nSvR,2),'color','k','LineWidth',4)
xlabel('E(condition) ./ E(ranked)')
ylabel('probability')
legend('Non-Shuffle Condition','Shuffle Condition')
axis tight
ylims = get(gca,'YLim');
% set(gca,'YLim',[0,0.15])
set(gca,'YLim',[0,ylims(2)+0.05])
for n = 1:size(s,2)
    scatter(mean(r{n}./e{n}),max(ylims)+0.02,200,'*','b')
    scatter(mean(r{n}./s{n}),max(ylims)+0.02,200,'*','k')
end

%%
oiE = {};
oiS = {};
oiR = {};
for n = 1:size(e,2)
    normfactor(n) = max([max(e{n}),max(s{n}),max(r{n})]);
    oiE{n} = e{n}./normfactor(n);
    oiS{n} = s{n}./normfactor(n);
    oiR{n} = r{n}./normfactor(n);
    
end


edges = [0:0.01:1];
nMotif = [];
nShuffle = [];
for n = 1:size(e,2)
    nE(:,n) = histc(oiE{n},edges);
    nS(:,n) = histc(oiS{n},edges);
    nR(:,n) = histc(oiR{n},edges);
end
figure;hold on
stairs(edges,sum(nR,2)./sum(sum(nR,1)),'color','r','LineWidth',4)
stairs(edges,sum(nS,2)./sum(sum(nS,1)),'color','b','LineWidth',4)
stairs(edges,sum(nE,2)./sum(sum(nE,1)),'color','k','LineWidth',4)
xlabel('Normalized Entropy')
ylabel('probability(Entropy)')
legend('Maximally Peaked','Shuffle Condition','Non-Shuffle Condition')
axis tight
ylims = get(gca,'YLim');
% set(gca,'YLim',[0,0.15])
set(gca,'YLim',[0,ylims(2)+0.05])
for n = 1:size(s,2)
    scatter(mean(oiE{n}),max(ylims)+0.02,200,'*','b')
    scatter(mean(oiS{n}),max(ylims)+0.02,200,'*','k')
    scatter(mean(oiR{n}),max(ylims)+0.02,200,'*','r')
end



edges = [0:0.01:1];
nMotif = [];
nShuffle = [];
for n = 1:size(e,2)
    nE(:,n) = histc(oiE{n} - oiR{n},edges);
    nS(:,n) = histc(oiS{n} - oiR{n},edges);
    %     nR(:,n) = histc(oiR{n},edges);
end
figure;hold on
% stairs(edges,sum(nR,2)./sum(sum(nR,1)),'color','r','LineWidth',4)
stairs(edges,sum(nS,2)./sum(sum(nS,1)),'color','b','LineWidth',4)
stairs(edges,sum(nE,2)./sum(sum(nE,1)),'color','k','LineWidth',4)
xlabel('Normalized Entropy (diff from minimum Entropy)')
ylabel('probability(Entropy)')
legend('Shuffle Condition','Non-Shuffle Condition')
axis tight
ylims = get(gca,'YLim');
% set(gca,'YLim',[0,0.15])
set(gca,'YLim',[0,ylims(2)+0.05])

for n = 1:size(s,2)
    scatter(mean(oiE{n}-oiR{n}),max(ylims)+0.02,200,'*','b')
    scatter(mean(oiS{n}-oiR{n}),max(ylims)+0.02,200,'*','k')
    %     scatter(mean(oiR{n}),max(ylims),200,'*','r')
end

%%


nstims = size(Meta_Vm,3);
nsubsample = size(Meta_Vm,1);
%%%%%%%%%%%%%%%%%%
% similarity among stims averaged across shuffles
ishuffle = 1;
istim = 1;
figure
hold on
plot(squeeze(Meta_Vm(ishuffle,:,istim)),'color','k')
plot(squeeze(Meta_VmHat_Ridge(ishuffle,istim,:)),'color','b')
plot(squeeze(Meta_VmHat_Lasso(ishuffle,istim,:)),'color','r')
axis tight
legend('black:Vm' , 'blue:ridge' , 'red:lasso')


%%%%%%%%%%%%%%%%%%%
%%%%% DEGREES OF FREEDOM AT CRITERIA UNDER EACH MODEL
%%%%%%%%%%%%%%%%%%%
DF_L = [];
DF_R = [];
for ishuffle = 1:nsubsample
    for istim = 1:nstims
        
        ialpha = 1;
        LassoFit = AllFits{ishuffle,ialpha,istim};
        lambdaind = max(find(LassoFit.MSE < Meta_MeanMSE(ishuffle,istim)));
        DF_L(ishuffle,istim) = LassoFit.DF(lambdaind);
        
        ialpha = 2;
        RidgeFit = AllFits{ishuffle,ialpha,istim};
        lambdaind = max(find(RidgeFit.MSE < Meta_MeanMSE(ishuffle,istim)));
        DF_R(ishuffle,istim) = RidgeFit.DF(lambdaind);
        
    end
end
figure;hold on
scatter(mean(DF_L,1),mean(DF_R,1),200,'k','fill')
axis tight
SetAxisUnity(gcf)
xlabel('DF at criteria under Lasso (mean across shuffle)')
ylabel('DF atcriteria under Ridge (mean across shuffle)')

%%%%%%%%%%%%%%%%%%%
%%%%% MEAN BETA COEFF ACROSS POPULATION AT EACH LAMBDA
%%%%% (avg across shuffles) then rank order and plot IMAGESC
%*******SHIT
% I can't do this if Lambda is different vector for each shuffle iteration
% maybe I can do some work with aligning L, but not on log scale gah
%%%%%%%%%%%%%%%%%%%

%- LASSO
ialpha = 1;
for ishuffle = 1:nsubsample
    for istim = 1:nstims
        LassoFit = AllFits{ishuffle,ialpha,istim};
        Lasso_Mat(ishuffle,istim,:) = LassoFit.Lambda;
        Lasso_Beta_Mat(ishuffle,istim,:,:) = squeeze(Meta_B(ishuffle,ialpha,istim,:,:));
        
    end
end

% - RIDGE
ialpha = 2;
for ishuffle = 1:nsubsample
    for istim = 1:nstims
        RidgeFit = AllFits{ishuffle,ialpha,istim};
        Ridge_Mat(ishuffle,istim,:) = RidgeFit.Lambda;
        Ridge_Beta_Mat(ishuffle,istim,:,:) = squeeze(Meta_B(ishuffle,ialpha,istim,:,:));
        
    end
end












figure;hold on
scatter(MCorr_Lasso_stim,MCorr_Ridge_stim,200,'k','fill')
SetAxisUnity(gcf)
xlabel('pairwise MCC under Lasso')
ylabel('pairwise MCC under Ridge')

MCorr_RvL = [];
for istim = 1:nstims
    R1 = matthewscorr(RidgeBCriteria(istim,:),LassoBCriteria(istim,:));
    R2 = matthewscorr(LassoBCriteria(istim,:),RidgeBCriteria(istim,:));
    MCorr_RvL(istim) = mean([R1,R2]);
    
end


edges = [-1:0.05:1];
n_RvL = histc(MCorr_RvL,edges);
n_Ridge = histc(MCorr_Ridge_stim,edges);
n_Lasso = histc(MCorr_Lasso_stim,edges);
figure; hold on
stairs(edges, n_RvL./size(MCorr_RvL,2), 'color','k','LineWidth',4)
stairs(edges, n_Ridge./size(MCorr_Ridge_stim,2), 'color','b','LineWidth',3)
stairs(edges, n_Lasso./size(MCorr_Lasso_stim,2), 'color','r','LineWidth',3)
legend('Ridge vs Lasso within stim','Ridge pairwise across stims','Lasso pairwise across stims')
xlabel('Mathews Correlation Coefficient')
ylabel('p(MCC)')

%%%%%%%%%%%%%%%%%%
% similarity across shuffles
nshuffles = 10;
istim = 1;

figure;hold on

for itrialshuffle = 1:nshuffles
    subplot(nshuffles/2,2,itrialshuffle)
    hold on
    plot(squeeze(Meta_Vm(itrialshuffle,:,istim)),'color','k')
    plot(sum(squeeze(Meta_VmHat_Ridge(itrialshuffle,istim,:)),2),'color','b')
    plot(sum(squeeze(Meta_VmHat_Lasso(itrialshuffle,istim,:)),2),'color','r')
    axis tight
    if i ==1
        legend('black:Vm' , 'blue:ridge' , 'red:lasso')
    end
end


DF_L = [];
DF_R = [];
for itrialshuffle = 1:nshuffles
    LassoFit = AllFits{itrialshuffle,1,istim};
    lambdaind = max(find(LassoFit.MSE < Meta_MeanMSE(itrialshuffle,istim)));
    DF_L(itrialshuffle) = LassoFit.DF(lambdaind);
    RidgeFit = AllFits{itrialshuffle,4,istim};
    lambdaind = max(find(RidgeFit.MSE < Meta_MeanMSE(itrialshuffle,istim)));
    DF_R(itrialshuffle) = RidgeFit.DF(lambdaind);
end
figure;hold on
scatter(DF_L,DF_R,200,'k','fill')
axis tight
SetAxisUnity(gcf)
xlabel('DF at criteria under Lasso')
ylabel('DF at criteria under Ridge')
title('each point is a different shuffle for one stimulus')

ialpha = 1;
istim = 1;
LassoFit = AllFits{itrialshuffle,ialpha,istim};
Lasso_B = squeeze(Meta_B(itrialshuffle,ialpha,istim,:,:));
lassoPlot(Lasso_B,LassoFit,'PlotType','CV');
lassoPlot(Lasso_B,LassoFit,'PlotType','Lambda');
set(gca,'YLim',[-1,1])

ialpha = 4;
RidgeFit = AllFits{itrialshuffle,ialpha,istim};
Ridge_B = squeeze(Meta_B(1,4,1,:,:));
lassoPlot(Ridge_B,RidgeFit,'PlotType','CV');
lassoPlot(Ridge_B,RidgeFit,'PlotType','Lambda');
set(gca,'YLim',[-1,1])



ialpha = 1; %for lasso
LassoBCriteria_analog = squeeze(Meta_Criteria_B(:,ialpha,istim,:));
LassoBCriteria = LassoBCriteria_analog;
LassoBCriteria(find(LassoBCriteria)) = 1;
MCorr_Lasso_stim = [];
Corrcoef_Lasso_stim = [];
pairind = 1;
for ishuffle1 = 1:nshuffles -1
    for ishuffle2 = ishuffle1 +1 :nshuffles
        R1 = matthewscorr(LassoBCriteria(ishuffle1,:), LassoBCriteria(ishuffle2,:));
        R2 = matthewscorr(LassoBCriteria(ishuffle2,:), LassoBCriteria(ishuffle1,:));
        MCorr_Lasso_stim(pairind) = mean([R1,R2]);
        
        
        c = corrcoef(LassoBCriteria_analog(ishuffle1,:), LassoBCriteria_analog(ishuffle2,:));
        Corrcoef_Lasso_stim(pairind) = c(idx);
        
        pairind = pairind+1;
    end
end
figure;imagesc(LassoBCriteria)
title('Lasso')
ylabel('shuffleID')
xlabel('unitID')
figure;
scatter(Corrcoef_Lasso_stim,MCorr_Lasso_stim,200,'k','fill')
xlabel('correlation coefficient')
ylabel('MCC')
title('under Lasso conditions pairwise shuffles')
SetAxisUnity(gcf)


idx = true(2);
idx = ~tril(idx);
ialpha = 4;
RidgeBCriteria_analog = squeeze(Meta_Criteria_B(:,ialpha,istim,:));
RidgeBCriteria = RidgeBCriteria_analog;
RidgeBCriteria(find(RidgeBCriteria)) = 1;
MCorr_Ridge_stim = [];
Corrcoef_Ridge_stim = [];
pairind = 1;
for ishuffle1 = 1:nshuffles -1
    for ishuffle2 = ishuffle1 +1 :nshuffles
        R1 = matthewscorr(RidgeBCriteria(ishuffle1,:), RidgeBCriteria(ishuffle2,:));
        R2 = matthewscorr(RidgeBCriteria(ishuffle2,:), RidgeBCriteria(ishuffle1,:));
        MCorr_Ridge_stim(pairind) = mean([R1,R2]);
        
        
        c = corrcoef(RidgeBCriteria_analog(ishuffle1,:), RidgeBCriteria_analog(ishuffle2,:));
        Corrcoef_Ridge_stim(pairind) = c(idx);
        
        pairind = pairind+1;
    end
end
figure;imagesc(RidgeBCriteria)
title('Ridge')
ylabel('shuffleID')
xlabel('unitID')
figure;
scatter(Corrcoef_Ridge_stim,MCorr_Ridge_stim,200,'k','fill')
xlabel('correlation coefficient')
ylabel('MCC')
title('under Ridge conditions pairwise shuffles')
SetAxisUnity(gcf)


figure;hold on
scatter(MCorr_Lasso_stim,MCorr_Ridge_stim,200,'k','fill')
SetAxisUnity(gcf)
xlabel('pairwise MCC under Lasso')
ylabel('pairwise MCC under Ridge')
title('across shuffles')

MCorr_RvL = [];
for itrialshuffle = 1:nshuffles
    R1 = matthewscorr(RidgeBCriteria(istim,itrialshuffle),LassoBCriteria(istim,itrialshuffle));
    R2 = matthewscorr(LassoBCriteria(istim,itrialshuffle),RidgeBCriteria(istim,itrialshuffle));
    MCorr_RvL(itrialshuffle) = mean([R1,R2]);
    
end


edges = [-1:0.05:1];
n_RvL = histc(MCorr_RvL,edges);
n_Ridge = histc(MCorr_Ridge_stim,edges);
n_Lasso = histc(MCorr_Lasso_stim,edges);
figure; hold on
stairs(edges, n_RvL./size(MCorr_RvL,2), 'color','k','LineWidth',4)
stairs(edges, n_Ridge./size(MCorr_Ridge_stim,2), 'color','b','LineWidth',3)
stairs(edges, n_Lasso./size(MCorr_Lasso_stim,2), 'color','r','LineWidth',3)
legend('Ridge vs Lasso within shuffle','Ridge pairwise across shuffles','Lasso pairwise across shuffles')
xlabel('Mathews Correlation Coefficient')
ylabel('p(MCC)')










VmHatdata = [];
for iterm = 1:size(B,1)
    VmHatdata(:,iterm) = B(iterm,bind)*thisPop(:,iterm);
end
VmHatdata;

Lambda_Vm_B_BINARY = MetaLambda_Vm_B;
Lambda_Vm_B_BINARY(find(Lambda_Vm_B_BINARY)) = 1;



MCorr_shuffle = [];
pairind = 1;
for istim1 = 1:nstims -1
    for istim2 = istim1 +1 :nstims
        R1 = matthewscorr(Lambda_Vm_B_BINARY(:,istim1), Lambda_Vm_B_BINARY(:,istim2));
        R2 = matthewscorr(Lambda_Vm_B_BINARY(:,istim2), Lambda_Vm_B_BINARY(:,istim1));
        MCorr_shuffle(pairind) = mean([R1,R2]);
        pairind = pairind+1;
    end
end

nstims = size(MetaLambda_Vm_B,3);
iterations = size(MetaLambda_Vm_B,1);
nclusters = size(MetaLambda_Vm_B,2);


SumDot_stim = [];
SumDot_shuffle = [];
for istim = 1:nstims
    
    thisstim = squeeze(MetaLambda_Vm_B(:,:,istim));
    
    pairind = 1;
    for ishuffle1 = 1:iterations-1
        for ishuffle2 = ishuffle1+1:iterations
            SumDot_shuffle(istim,pairind) = sum(thisstim(ishuffle1,:).*thisstim(ishuffle2,:));
            pairind = pairind+1;
        end
    end
end
for ishuffle = 1:iterations
    
    thisshuffle = squeeze(MetaLambda_Vm_B(ishuffle,:,:));
    
    pairind = 1;
    for istim1 = 1:nstims-1
        for istim2 = istim1+1:nstims
            SumDot_stim(ishuffle,pairind) = sum(thisshuffle(:,istim1).*thisshuffle(:,istim2));
            pairind = pairind+1;
        end
    end
    
end

estimate_overlap_stims = mean(SumDot_stim,2)
estimate_overlap_shuffles = mean(SumDot_shuffle,2)

[h,p] = kstest2(estimate_overlap_stims,estimate_overlap_shuffles)

%%
lassoPlot(B,FitInfo,'PlotType','CV');
istim = 1;

hfig_scatt = figure;
hold on
hfig_predict = figure;
hold on

colors = {'r','k','m','b'};
lassoinds = [99,max(find(FitInfo.MSE < MeanMSE(istim))),FitInfo.Index1SE,2];
for ind = 1:size(lassoinds,2)
    
    
    VmHatdata = [];
    for iterm = 1:size(B,1)
        VmHatdata(:,iterm) = B(iterm,lassoinds(ind))*thisPop(:,iterm);
    end
    
    figure(hfig_scatt);
    hold on
    [x,i] = sort(B(:, lassoinds(ind)));
    scatter([1:length(x)],x,colors{ind})
    
    figure(hfig_predict);
    hold on
    plot(sum(VmHatdata,2),'color',colors{ind})
end


figure(hfig_scatt);
set(gca,'YLim',[-1,1])
legend('max lambda','criteria lambda','1SE of min lambda','min lambda')
figure(hfig_predict);
legend('max lambda','criteria lambda','1SE of min lambda','min lambda')

%% after getting estimate overlap for all samples... and copy/paste from excel doc

% tmp = {};
% for n = 1:size(ShuffleOverlap,2)
%     tmp{n} = cell2mat(ShuffleOverlap(:,n));
% end
% ShuffleOverlap = tmp;

m = [];
s = {};
nPredictors = [30,56,56,41,41,41,41,41];
for n = 1:size(ShuffleOverlap,2)
    m(:,n) = (MotifOverlap(:,n)./nPredictors(n));
    s{n} = (ShuffleOverlap{n}./nPredictors(n));
end

for n = 1:size(s,2)
    [h,p] = kstest2(s{n},m(:,n));
    isDifferent(n) = h;
    howDifferent(n) = p;
    
end

edges = [0:.01:1];
nMotif = [];
nShuffle = [];
for n = 1:size(s,2)
    nShuffle(:,n) = histc(s{n},edges);
    nMotif(:,n) = histc(m(:,n),edges);
end
figure;hold on
stairs(edges,sum(nMotif,2)./sum(sum(nMotif,1)),'color','k','LineWidth',4)
stairs(edges,sum(nShuffle,2)./sum(sum(nShuffle,1)),'color','b','LineWidth',4)
xlabel('Overlap Index')
ylabel('probability(OI)')
legend('Motif Condition','Shuffle Condition')
set(gca,'YLim',[0,0.18])
ylims = get(gca,'YLim');
for n = 1:size(s,2)
    scatter(mean(s{n}),max(ylims),200,'*','b')
    scatter(mean(m(:,n)),max(ylims),200,'*','k')
end
set(gca,'YLim',[0,0.25])


alls = [];
allm = [];
for n = 1:size(s,2)
    alls = [alls;s{n}];
    allm = [allm;m(:,n)];
end
meanMotifOverlap = mean(allm)
getCDFconf(allm,95)
meanShuffleOverlap = mean(alls)
getCDFconf(alls,95)

oiNormShuffle = {};
oiNormMotif = [];
for n = 1:size(ShuffleOverlap,2)
    normfactor(n) = max([max(ShuffleOverlap{n}),max(MotifOverlap(:,n))]);
    oiNormShuffle{n} = ShuffleOverlap{n}./normfactor(n);
    oiNormMotif(:,n) = MotifOverlap(:,n)./normfactor(n);
    
end


edges = [0:0.01:1];
nMotif = [];
nShuffle = [];
for n = 1:size(ShuffleOverlap,2)
    nShuffle(:,n) = histc(oiNormShuffle{n},edges);
    nMotif(:,n) = histc(oiNormMotif(:,n),edges);
end
figure;hold on
stairs(edges,sum(nMotif,2)./sum(sum(nMotif,1)),'color','k','LineWidth',4)
stairs(edges,sum(nShuffle,2)./sum(sum(nShuffle,1)),'color','b','LineWidth',4)
xlabel('Normalized Overlap Index')
ylabel('probability(OI)')
legend('Motif Condition','Shuffle Condition')
set(gca,'YLim',[0,0.18])
ylims = get(gca,'YLim');
for n = 1:size(s,2)
    scatter(mean(oiNormShuffle{n}),max(ylims),200,'*','b')
    scatter(mean(oiNormMotif(:,n)),max(ylims),200,'*','k')
end
set(gca,'YLim',[0,0.21])

%% get spectral* entropy of spatial distribution beta weights
BirdID = {'B1075'
    'B1083'
    'B1083'
    'B1087'
    'B1087'
    'B1112'
    'B1112'
    'B1112'};

Sample = {'ClusterID_Overlap_MetaVectors_CellID1.mat'
    'ClusterID_Overlap_MetaVectors_CellID1.mat'
    'ClusterID_Overlap_MetaVectors_CellID2.mat'
    'ClusterID_Overlap_MetaVectors_CellID1.mat'
    'ClusterID_Overlap_MetaVectors_CellID2.mat'
    'ClusterID_Overlap_MetaVectors_CellID1.mat'
    'ClusterID_Overlap_MetaVectors_CellID2.mat'
    'ClusterID_Overlap_MetaVectors_CellID3.mat'};
%%
for pairID = 1:size(BirdID,1);
    cd(['/Users/kperks/mnt/cube/Ice/kperks/' BirdID{pairID} '/MetaMat/'])
    load(Sample{pairID},'MetaClusterIDMat','MetaBetaWeightMatrix')
    iterations = size(MetaBetaWeightMatrix,1);
    nstims = size(MetaBetaWeightMatrix,3);
    
    load('MVRegressPrep.mat','include_inds','isNarrow')
    load('metatoes.mat','allchanxz','allclusterid')
    
    these_locations = allchanxz(include_inds);
    locationVec = unique(these_locations);
    BetaVec = zeros(iterations,nstims,size(locationVec,1),1);
    BetaVec_shuffle = zeros(iterations,nstims,size(locationVec,1),1);
    Entropy = nan(iterations,nstims);
    Entropy_shuffle = nan(iterations,nstims);
    Entropy_ranked = nan(iterations,nstims);
    for i = 1:iterations
        
        for istim = 1:nstims
            stimB = abs(squeeze(MetaBetaWeightMatrix(i,:,istim)));
            
            these_locations = allchanxz(include_inds);
            
            shuffleinds = randperm(size(stimB,2));
            stimB_shuffle = stimB(shuffleinds);
            %         [x,i] = sort(stimB);
            nloc = [];
            for l = 1:size(locationVec,1)
                linds = find(these_locations == locationVec(l));
                nloc(l) = size(linds,1);
                if ~isempty(linds)
                    BetaVec(i,istim,l) = mean(stimB(linds));
                    BetaVec_shuffle(i,istim,l) = mean(stimB_shuffle(linds));
                end
            end
            
            %%%%%Shuffle control
            d = squeeze(BetaVec_shuffle(i,istim,:));
            d=d/sum(d); %Normalization
            %Entropy Calculation
            logd = log2(d);
            logd(isinf(logd)) = 0;
            Entropy_shuffle(i,istim) = -sum(d.*logd)/log2(size(d,1));
            
            %%%%%Non-shuffle condition
            d = squeeze(BetaVec(i,istim,:));
            d=d/sum(d); %Normalization
            %Entropy Calculation
            logd = log2(d);
            logd(isinf(logd)) = 0;
            Entropy(i,istim) = -sum(d.*logd)/log2(size(d,1));
            
            %%%%%Rank Order Beta by location to get max peakiness
            sortlocN = nloc; %sort(nloc,'descend'); %
            sortB = sort(stimB,'descend');
            distB = [];
            for l = 1:size(locationVec,1)
                if ~isempty(find(these_locations == locationVec(l)))
                    distB(l) = mean(sortB(1:sortlocN(l)));
                    sortB = sortB(sortlocN(l)+1:end);
                end
            end
            d_rank = distB ./ sum(distB);%Normalization
            %Entropy Calculation
            logd = log2(d_rank);
            logd(isinf(logd)) = 0;
            Entropy_ranked(i,istim) = -sum(d_rank.*logd)/log2(size(d_rank,2));
            
        end
    end
    
    Entropy_meta{pairID} = mean(Entropy,1)';
    Entropy_meta_Shuffle{pairID} = mean(Entropy_shuffle,1)';
    Entropy_meta_Ranked{pairID} = mean(Entropy_ranked,1)';
end

figure;plot(locationVec,d)
figure;plot(locationVec,d_rank)
% after getting estimate entropy for all samples... and copy/paste from excel doc
%for each sample, plot relative to minimum possible entropy(most peaky possible)
e = Entropy_meta;
s = Entropy_meta_Shuffle;
r = Entropy_meta_Ranked;

alls = [];
alle = [];
allr = [];
for n = 1:size(s,2)
    alls = [alls;s{n}];
    alle = [alle;e{n}];
    allr = [allr;r{n}];
end
meanEntropy = mean(alle)
getCDFconf(alle,95)
meanEntropyShuffle = mean(alls)
getCDFconf(alls,95)
meanEntropyRanked = mean(allr)
getCDFconf(allr,95)
meanEntropy = mean(alle-allr)
getCDFconf(alle-allr,95)
meanEntropyShuffle = mean(alls-allr)
getCDFconf(alls-allr,95)

oiE = {};
oiS = {};
oiR = {};
for n = 1:size(e,2)
    normfactor(n) = max([max(e{n}),max(s{n}),max(r{n})]);
    oiE{n} = e{n}./normfactor(n);
    oiS{n} = s{n}./normfactor(n);
    oiR{n} = r{n}./normfactor(n);
    
end


edges = [0:0.01:1];
nMotif = [];
nShuffle = [];
for n = 1:size(e,2)
    nE(:,n) = histc(oiE{n},edges);
    nS(:,n) = histc(oiS{n},edges);
    nR(:,n) = histc(oiR{n},edges);
end
figure;hold on
stairs(edges,sum(nR,2)./sum(sum(nR,1)),'color','r','LineWidth',4)
stairs(edges,sum(nS,2)./sum(sum(nS,1)),'color','b','LineWidth',4)
stairs(edges,sum(nE,2)./sum(sum(nE,1)),'color','k','LineWidth',4)
xlabel('Normalized Entropy')
ylabel('probability(Entropy)')
legend('Maximally Peaked','Shuffle Condition','Non-Shuffle Condition')
axis tight
ylims = get(gca,'YLim');
% set(gca,'YLim',[0,0.15])
set(gca,'YLim',[0,ylims(2)+0.05])
for n = 1:size(s,2)
    scatter(mean(oiE{n}),max(ylims)+0.02,200,'*','b')
    scatter(mean(oiS{n}),max(ylims)+0.02,200,'*','k')
    scatter(mean(oiR{n}),max(ylims)+0.02,200,'*','r')
end



edges = [0:0.01:1];
nMotif = [];
nShuffle = [];
for n = 1:size(e,2)
    nE(:,n) = histc(oiE{n} - oiR{n},edges);
    nS(:,n) = histc(oiS{n} - oiR{n},edges);
    %     nR(:,n) = histc(oiR{n},edges);
end
figure;hold on
% stairs(edges,sum(nR,2)./sum(sum(nR,1)),'color','r','LineWidth',4)
stairs(edges,sum(nS,2)./sum(sum(nS,1)),'color','b','LineWidth',4)
stairs(edges,sum(nE,2)./sum(sum(nE,1)),'color','k','LineWidth',4)
xlabel('Normalized Entropy (diff from minimum Entropy)')
ylabel('probability(Entropy)')
legend('Shuffle Condition','Non-Shuffle Condition')
axis tight
ylims = get(gca,'YLim');
% set(gca,'YLim',[0,0.15])
set(gca,'YLim',[0,ylims(2)+0.05])

for n = 1:size(s,2)
    scatter(mean(oiE{n}-oiR{n}),max(ylims)+0.02,200,'*','b')
    scatter(mean(oiS{n}-oiR{n}),max(ylims)+0.02,200,'*','k')
    %     scatter(mean(oiR{n}),max(ylims),200,'*','r')
end
%% Wide / Narrow mean population level activity of only predictive untis
for istim = 1:nstims
    stimB = BetaWeightMatrix(:,istim);
    [x,i] = sort(stimB);
    
    %     increment = signifCorr_Increment(istim);
    increment = size(corrInc,2);
    
    stimVmHat = squeeze(VmHat(istim,:,:));
    stimVm = Vm(:,istim);
    stimPVm = squeeze(PVm(istim,:,:));
    
    
    Bneg = i(1:increment);
    %         Bneg = [];
    Bpos = i(end-increment+1:end);
    Binds = [Bneg;Bpos];
    thisB = stimB(Binds,1);
    
    thisPopHat_Narrow = sum(stimVmHat(:,intersect(Binds,find(NarrowInds ==1))),2);
    thisPopHat_Wide = sum(stimVmHat(:,intersect(Binds,find(NarrowInds ==0))),2);
    
    thisPop_Narrow = sum(stimPVm(:,intersect(Binds,find(NarrowInds ==1))),2);
    thisPop_Wide = sum(stimPVm(:,intersect(Binds,find(NarrowInds ==0))),2);
    
    thisPopHat_Bneg = sum(stimVmHat(:,Bneg),2);
    thisPopHat_Bpos = sum(stimVmHat(:,Bpos),2);
    
    thisPop_Bneg = sum(stimPVm(:,Bneg),2);
    thisPop_Bpos = sum(stimPVm(:,Bpos),2);
    
    %         thisPop = sum(stimVmHat(:,i(end-increment+1:end)),2);
    %         thisPop = (thisPop-min(thisPop));
    %         thisPop = thisPop./max(thisPop);
    %         popInc(istim,:,increment) = thisPop;
    %         varInc(istim,increment) = mean(var([thisPop';stimVm']));
    %         c = corrcoef([thisPop';stimVm']');
    %         corrInc(istim,increment) =c(idx);
    %         sumR2Inc(istim,increment) = sum(power(diff([thisPop';stimVm']),2));
    
    
end

figure;hold on
xtime = ([1:size(thisPopHat_Narrow,1)]*dt)+calcstart*dt;
line(xtime,sum([thisPopHat_Wide,thisPopHat_Narrow],2),'color','k')
line(xtime,thisPopHat_Narrow,'color','r')
line(xtime,thisPopHat_Wide,'color','b')
line(xtime,stimVm,'color','m');
axis tight
legend('sum wide plus narrow (B=trained)','Narrow Unit VmHat Sum','Wide Unit VmHat Sum','Vm')
% legend('Narrow Unit VmHat Sum','Wide Unit VmHat Sum','Vm')
title(['stimulus ' num2str(istim) ' B1083 n=' num2str(increment*2) ' units in population was necessary'])

figure;hold on
xtime = ([1:size(thisPop_Narrow,1)]*dt)+calcstart*dt;
line(xtime,sum([thisPop_Wide,thisPop_Narrow],2),'color','k')
line(xtime,thisPop_Narrow,'color','r')
line(xtime,thisPop_Wide,'color','b')
line(xtime,stimVm,'color','m');
axis tight
% legend('sum wide plus narrow (B=1)','Narrow Unit B=1 Sum','Wide Unit B=1 Sum','Vm')
legend('Narrow Unit B=1 Sum','Wide Unit B=1 Sum','Vm')
title(['stimulus ' num2str(istim) ' B1083 n=' num2str(increment*2) ' units in population was necessary'])

figure;hold on
xtime = ([1:size(thisPopHat_Bneg,1)]*dt)+calcstart*dt;
line(xtime,sum([thisPopHat_Bpos,thisPopHat_Bneg],2),'color','k')
line(xtime,thisPopHat_Bneg,'color','r')
line(xtime,thisPopHat_Bpos,'color','b')
line(xtime,stimVm,'color','m');
axis tight
legend('sum wide plus narrow (B=trained)','Beta Neg Unit VmHat Sum','Beta Wide Unit VmHat Sum','Vm')
% legend('Beta Neg Unit VmHat Sum','Beta Wide Unit VmHat Sum','Vm')
title(['stimulus ' num2str(istim) ' B1083 n=' num2str(increment*2) ' units in population was necessary'])

figure;hold on
xtime = ([1:size(thisPop_Bneg,1)]*dt)+calcstart*dt;
% line(xtime,sum([thisPop_Bpos,thisPop_Bneg],2),'color','k')
line(xtime,thisPop_Bneg,'color','r')
line(xtime,thisPop_Bpos,'color','b')
line(xtime,stimVm,'color','m');
axis tight
% legend('sum wide plus narrow (B=1)','Beta Neg Unit B=1 Sum','Beta Wide Unit B=1 Sum','Vm')
legend('Beta Neg Unit B=1 Sum','Beta Wide Unit B=1 Sum','Vm')
title(['stimulus ' num2str(istim) ' B1087 n=' num2str(increment*2) ' units in population '])

echo_corrBnegBpos = corrcoef(thisPopHat_Bneg,thisPopHat_Bpos)
echo_corrBnegVm = corrcoef(thisPopHat_Bneg,stimVm)
echo_corrBposVm = corrcoef(thisPopHat_Bpos,stimVm)
echo_corrSumPopVm = corrcoef(sum([thisPopHat_Bpos,thisPopHat_Bneg],2),stimVm)

%%

location_mat = zeros(size(unique(loc),1),nstims);
locationNeg_mat = zeros(size(unique(loc),1),nstims);
locationPos_mat = zeros(size(unique(loc),1),nstims);
loc_array = unique(loc);
for istim = 1:nstims
    stimB = BetaWeightMatrix(:,istim);
    [x,i] = sort(stimB);
    
    %     this_increment = round(cell_signifCorr(istim)/100*nclusters)/2;
    %     increment = find(percentTerms == cell_signifCorr(istim));
    increment = signifCorr_Increment(istim);
    if ~isnan(increment)
        Bneg = i(1:increment);
        %         Bneg = [];
        Bpos = i(end-increment+1:end);
        Binds = [Bneg;Bpos];
        %         these_inds = i([i(1:increment),i(end-increment+1:end)],1);
        these_locations = loc(Binds);
        %         these_locations_neg = loc(Bneg);
        %         these_locations_pos = loc(Bpos);
        these_locations_neg = loc(intersect(Binds,find(NarrowInds==1)));
        these_locations_pos = loc(intersect(Binds,find(NarrowInds==0)));
        for ib = 1:size(Binds,1)
            location_mat(find(loc_array == these_locations(ib)),istim) = 1;
        end
        for ib = 1:size(these_locations_neg,1)
            locationNeg_mat(find(loc_array == these_locations_neg(ib)),istim) = 1;
        end
        for ib = 1:size(these_locations_pos,1)
            locationPos_mat(find(loc_array == these_locations_pos(ib)),istim) = 1;
        end
        
    end
end
figure;imagesc(location_mat)
set(gca,'YTick',[1:4:size(loc_array,1)],'YTickLabel',loc_array([1:4:size(loc_array,1)]))
set(gca,'XTick',[1:nstims],'XTickLabel',metatoes_songlist)
ylabel('location along penetration')
xlabel('PseudoSong label')
title('significance by Correlation')
colorbar
set(gcf,'Position',[17   233   299   859])

figure;imagesc(locationNeg_mat)
set(gca,'YTick',[1:4:size(loc_array,1)],'YTickLabel',loc_array([1:4:size(loc_array,1)]))
set(gca,'XTick',[1:nstims],'XTickLabel',metatoes_songlist)
ylabel('location along penetration')
xlabel('PseudoSong label')
% title('(-) Beta (signif by Corr)')
title('Narrow Units (signif by Corr)')
colorbar
set(gcf,'Position',[17   233   299   859])

figure;imagesc(locationPos_mat)
set(gca,'YTick',[1:4:size(loc_array,1)],'YTickLabel',loc_array([1:4:size(loc_array,1)]))
set(gca,'XTick',[1:nstims],'XTickLabel',metatoes_songlist)
ylabel('location along penetration')
xlabel('PseudoSong label')
title('Wide Units (signif by Corr)')
colorbar
set(gcf,'Position',[17   233   299   859])

figure;hold on
for istim = 1:nstims
    subplot(nstims/2,2,istim)
    ylabel(['corr stim ' num2str(istim)])
    line((([1:size(corrInc,2)]*2)/nclusters)*100,corrInc(istim,:))
    xlims = get(gca,'XLim');
    line(xlims,[mean(pairwiseTrialCorr{istim}),mean(pairwiseTrialCorr{istim})],'color','r')
    set(gca,'YLim',[0,1])
end
set(gcf,'Position',[1361         298         491         777])


%% Meta Analysis
% metavars here saved as MetaData_PredictionAnalysis.mat in mnt/cube/Ice/kperks

a = {};
%copy and paste data into a from excel and then process with:

a = a(~cellfun(@isempty,a));
a = a(~cellfun(@ischar,a));

MeanTrialCorr = cell2mat(a);
MaxCorrPredict = cell2mat(a);
PercentSignifTerms = cell2mat(a);

MaxCorrPredict_Shuffled = cell2mat(a);
PercentSignifTerms_Shuffled = cell2mat(a);

median(MeanTrialCorr)
getCDFconf(MeanTrialCorr,95)

median(MaxCorrPredict)
getCDFconf(MaxCorrPredict,95)

median(MaxCorrPredict_Shuffled)
getCDFconf(MaxCorrPredict_Shuffled,95)

edges = [-1:0.1:1];
[n1,bins] = histc(MeanTrialCorr,edges);
[n2,bins] = histc(MaxCorrPredict,edges);
[n3,bins] = histc(MaxCorrPredict_Shuffled,edges);

figure;hold on
stairs(edges,n1,'color','k')
stairs(edges,n2,'color','b')
stairs(edges,n3,'color','r')

mean(MaxCorrPredict_Shuffled./MaxCorrPredict)
getCDFconf((MaxCorrPredict_Shuffled./MaxCorrPredict),95)

figure;hold on
scatter(MaxCorrPredict,MaxCorrPredict_Shuffled,200,'k','fill')
set(gca,'YLim',[0.5,1],'XLim',[0.5,1])
ylabel('Maximum (Shuffled) Prediction Correlation')
xlabel('Maximum Prediction Correlation')

mean(PercentSignifTerms_Shuffled./PercentSignifTerms)
getCDFconf((PercentSignifTerms_Shuffled./PercentSignifTerms),95)

figure;hold on
scatter(PercentSignifTerms,PercentSignifTerms_Shuffled,200,'k','fill')
set(gca,'YLim',[0,100],'XLim',[0,100])
ylabel('Percent (Shuffled) Prediction Terms')
xlabel('Percent Prediction Terms')

median(PercentSignifTerms_Shuffled)
getCDFconf(PercentSignifTerms_Shuffled,95)

median(PercentSignifTerms)
getCDFconf(PercentSignifTerms,95)

edges = [0:5:100];
[n1,bins] = histc(PercentSignifTerms,edges);
[n2,bins] = histc(PercentSignifTerms_Shuffled,edges);
figure;hold on
stairs(edges,n1,'color','k')
stairs(edges,n2,'color','r')

