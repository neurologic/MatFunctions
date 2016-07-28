
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
%
Lalpha = 1;
Ralpha = 2;

DF_Ratio = [];
MCorr_lasso_stim_ALL = [];
MCorr_lasso_subsample_ALL = [];
MCorr_ridge_stim_ALL = [];
MCorr_ridge_subsample_ALL = [];
    
for pairID = 1:size(BirdID,1);
    fileID = fopen('matlog.txt','w');
    t = datestr(datetime('now'));
    fprintf(fileID,'%s %i %s\n','on pair ID', pairID, t);
    fclose(fileID);
%     cd(['/Users/kperks/mnt/cube/Ice/kperks/' BirdID{pairID} '/MetaMat/'])
    cd(['/mnt/cube/Ice/kperks/' BirdID{pairID} '/MetaMat/'])
    
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
