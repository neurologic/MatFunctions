cd('/mnt/cube/Ice/kperks/B1112/MetaMat/')
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

CellID = 3; %which intracellular cell

Pop = SignalFiltR;
calcstart = size(Pop,3)/6;
nstims = size(SignalFiltR,2);
nclusters = size(SignalFiltR,1);

stiminds = [1:nstims];

subset_win = [calcstart+1:(2*calcstart)];
% subset_win = [calcstart+1:size(Pop,3)];

nsamps = size(subset_win,2);

Vm = [];

  pairwiseTrialCorr = [];
    pairwiseTrialXCorr = [];
  
        
    idx = true(2);
    idx = ~tril(idx);
    
motifind = 1;
for istim = 1:nstims
    
    for imot = 1:5
        subset_win = [(imot*calcstart)+1:((imot+1)*calcstart)];
        
        
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
        
        Vmdata = mean(IntracellularData{CellID}.stim_data{istim}(:,1:min(trialsamps)),1);
        Vmdata = Vmdata(subset_win);
        Vmdata = Vmdata - globaloffset;
        Vmdata = Vmdata ./ globalnormfactor;
        %     Vm(:,istim) = fliplr(Vmdata)';
        Vm(:,motifind) = Vmdata';
        
        thisVm = IntracellularData{CellID}.stim_data{istim}(:,1:min(trialsamps));
        thisVm = thisVm - repmat(globaloffset,size(thisVm,1),size(thisVm,2));
        thisVm = thisVm ./ repmat(globalnormfactor,size(thisVm,1),size(thisVm,2));
        thisVm = thisVm(:,subset_win);
        %correlation metric for goodness of match to Vm
        pairind = 1;
        for itrial1 = 1:size(thisVm,1)-1
            for itrial2 = itrial1+1:size(thisVm,1)
                
                c = corrcoef(thisVm([itrial1,itrial2],:)');
                pairwiseTrialCorr{motifind}(pairind) = c(idx);
                
                [c,lags] = xcorr(thisVm(itrial1,:),thisVm(itrial2,:),'coeff');
                pairwiseTrialXCorr{motifind}(pairind) = max(c);
                
                pairind = pairind+1;
            end
        end
        
        motifind = motifind+1;
    end
end
MeanCorrs = cellfun(@(x) mean(x),pairwiseTrialCorr)';
MeanXCorrs = cellfun(@(x) mean(x),pairwiseTrialXCorr)';
%
nmotifs = size(Vm,2);

 MetaClusterIDMat = [];
 MetaBetaWeightMatrix = [];
for trialsuffleid = 1:10;
    % load('MVRegressPrep_Shuffled.mat','SignalFiltR','include_inds','isNarrow','ExtracellStimNames','IntracellStimNames')
    %     load('MVRegressPrep.mat','SignalFiltR','include_inds','isNarrow','ExtracellStimNames','IntracellStimNames')
    load(['MVRegressPrep_Subsample10Trial_' num2str(trialsuffleid) '.mat'],'SignalFiltR','include_inds','isNarrow')
    
    fileID = fopen('matlog.txt','w');
        t = datestr(datetime('now'));
        fprintf(fileID,'%s %i %s\n','on shuffle rep', trialsuffleid, t);
        fclose(fileID);
    Pop = SignalFiltR;
    
    %     Vm = zeros(nsamps,nstims);
    BetaWeightMatrix = zeros(nclusters,nmotifs);
    PVm = zeros(nmotifs,nsamps,nclusters);
    VmHat = zeros(nmotifs,nsamps,nclusters);
    SignalPower = zeros(nsamps,nmotifs);
    
    % get global normalization factor for each cluster across stims
    % (already offset to 0)
    tmpmax = abs((max(max(Pop,[],3),[],2)));
    tmpmin = abs((min(min(Pop,[],3),[],2)));
    cluster_globalnorm = max([tmpmax,tmpmin],[],2);
    
    motifind = 1;
    
    for istim = 1:nstims
                [wav,fs] = audioread(['/mnt/cube/Ice/6sPseudoSongs/' IntracellStimNames{istim} '.wav']);

%         [wav,fs] = audioread(['/Users/kperks/mnt/cube/Ice/6sPseudoSongs/' IntracellStimNames{istim} '.wav']);
        y = hilbert(wav);
        env = abs(y);
        dnsampV = round([1:size(env,1)/size(Pop,3):size(env,1)]);
        dnsampWav = env(dnsampV,:);
        
        for imot = 1:5
            subset_win = [(imot*calcstart)+1:((imot+1)*calcstart)];
            
            Pop_stim = squeeze(Pop(:,stiminds(istim),subset_win));
            Pop_stim = Pop_stim';
            
            PVmdata = [];
            for iterm = 1:size(Pop_stim,2)
                normterm = cluster_globalnorm(iterm) ;
                PVmdata(:,iterm) = Pop_stim(:,iterm)./normterm;
            end
            PVmdata(isnan(PVmdata)) = 0;
            PVm(motifind,:,:) = PVmdata;
            
            SignalPower(:,motifind) = dnsampWav(subset_win,1);
            
            motifind = motifind+1;
        end
    end
    
  
    preTrainingPopCorr = [];
    preTrainingPopXCorr = [];
    
    for istim = 1:nmotifs
        
        thisVm = Vm(:,istim);
        thisPop = squeeze(PVm(istim,:,:));
        
        
        c = corrcoef(thisVm',sum(thisPop,2)');
        preTrainingPopCorr{istim} = c(idx);
        c = xcorr(thisVm',sum(thisPop,2)','coeff');
        preTrainingPopXCorr{istim} = max(c);
        
        [B] = regress(thisVm,thisPop);
        
        BetaWeightMatrix(:,istim) = B;
        
        VmHatdata = [];
        for iterm = 1:size(B,1)
            VmHatdata(:,iterm) = B(iterm)*thisPop(:,iterm);
        end
        VmHat(istim,:,:) = VmHatdata;
        
        %         figure;hold on
        %         xtime = [1:size(VmHatdata,1)]*dt;
        %         line(xtime,sum(VmHatdata,2)')
        %         line(xtime,thisVm,'color','r')
        % %         legend('VmHat trained only on intro motif','Vm')
        %         set(gcf,'Position',[298         639        1306         459])
        %         title(['songID ' num2str(istim)])
        
    end
    
    %
    popInc = [];
    corrInc = [];
    SignalVmCorr = [];
    SignalPVmCorr = [];
    for istim = 1:nmotifs
        stimB = BetaWeightMatrix(:,istim);
        [x,i] = sort(stimB);
        
        stimVmHat = squeeze(VmHat(istim,:,:));
        stimVm = Vm(:,istim);
        stimPVm = squeeze(PVm(istim,:,:));
        
        stimPower = SignalPower(:,istim);
        c = corrcoef([stimVm';stimPower']');
        SignalVmCorr(istim) = c(idx);
        c = corrcoef([sum(stimPVm,2)';stimPower']');
        SignalPVmCorr(istim) = c(idx);
        
        for increment = 1:floor(size(stimB,1)/2)
            Bneg = i(1:increment);
            Bpos = i(end-increment+1:end);
            Binds = [Bneg;Bpos];
            thisB = stimB(Binds,1);
            thisPop = sum(stimVmHat(:,Binds),2);
            popInc(istim,:,increment) = thisPop;
            c = corrcoef([thisPop';stimVm']');
            corrInc(istim,increment) =c(idx);
        end
        
    end
    max_corr = max(corrInc')';
    echo_SignalVmCorr = SignalVmCorr';
    echo_SignalPVmCorr = SignalPVmCorr';
    %
    
    %     figure;hold on
    %     for istim = 1:nstims
    %         subplot(nstims/2,2,istim)
    %         ylabel(['corr stim ' num2str(istim)])
    %         line((([1:size(corrInc,2)]*2)/nclusters)*100,corrInc(istim,:))
    %         xlims = get(gca,'XLim');
    %         line(xlims,[mean(pairwiseTrialCorr{istim}),mean(pairwiseTrialCorr{istim})],'color','r')
    %         set(gca,'YLim',[0,1])
    %     end
    %     set(gcf,'Position',[1361         298         491         777])
    %
    
    percentTerms = (([1:size(corrInc,2)]*2)/nclusters)*100;
    nIncrement = [1:size(corrInc,2)];
    
    cell_signifCorr = nan(1,nstims);
    signifCorr_Increment = nan(1,nstims);
    
    

FractionPredictors = 0.15;
nPredictors = floor(FractionPredictors *nclusters);
    ClusterIDMat = [];
    for istim = 1:nmotifs
        stimB = BetaWeightMatrix(:,istim);
        [x,i] = sort(stimB);
        
        set_increment = nPredictors;
        %         set_increment = signifCorr_Increment(istim);
        
        Bneg = i(1:set_increment);
        %         Bneg = [];
        Bpos = i(end-set_increment+1:end);
        Binds = [Bneg;Bpos];
        %         these_inds = i([i(1:increment),i(end-increment+1:end)],1);
        
        clusterIDvector = zeros(nclusters,1);
        clusterIDvector(Binds) = 1;
        ClusterIDMat = [ClusterIDMat,clusterIDvector];
        
        if ~isempty(find(corrInc(istim,:)>mean(pairwiseTrialCorr{istim})))
            thresholdInd = min(find(corrInc(istim,:)>mean(pairwiseTrialCorr{istim})));
            cell_signifCorr(istim) = percentTerms(thresholdInd);
            signifCorr_Increment(istim) = nIncrement(thresholdInd);
            
        else
            cell_signifCorr(istim) = nclusters/nclusters*100;
            signifCorr_Increment(istim) = floor(nclusters/2);
            %                 clusterIDvector = ones(nclusters,1);
            %                 ClusterIDMat = [ClusterIDMat,clusterIDvector];
        end
        
    end
    
    
    
    echo_percentSignif = cell_signifCorr';
    echo_nTerms = signifCorr_Increment'*2;
    
    
    MetaClusterIDMat(trialsuffleid,:,:) = ClusterIDMat;
    MetaBetaWeightMatrix(trialsuffleid,:,:) = BetaWeightMatrix;
    
    % figure;plot(sum(squeeze(PVm(1,:,:)),2)')
    % title('subset 10 trials for each unit response')
    % title('subset 10 trials for each unit response before sum PVm stim1')
    % xlabel('unit ID')
    % ylabel('sum activity stim 1')
end


save(['ClusterID_Overlap_MetaVectors_CellID' num2str(CellID) '.mat'],'MetaClusterIDMat','MetaBetaWeightMatrix')
