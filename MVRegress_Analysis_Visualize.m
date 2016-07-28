%%
load('ClusterID_Overlap_MetaVectors_CellID1.mat','MetaClusterIDMat','MetaBetaWeightMatrix')

nstims = size(MetaClusterIDMat,3);
iterations = size(MetaClusterIDMat,1);
nclusters = size(MetaClusterIDMat,2);

FractionPredictors = 0.15;
nPredictors = floor(FractionPredictors *nclusters);

for ishuffle = 1:iterations
    ClusterIDMat = [];
    for istim = 1:nstims
        stimB = squeeze(MetaBetaWeightMatrix(ishuffle,:,istim)); %BetaWeightMatrix(:,istim);
        [x,ishuffle] = sort(stimB);
        
        set_increment = nPredictors;
        %         set_increment = signifCorr_Increment(istim);
        
        Bneg = ishuffle(1:set_increment);
        %         Bneg = [];
        Bpos = ishuffle(end-set_increment+1:end);
        Binds = [Bneg;Bpos];
        %         these_inds = i([i(1:increment),i(end-increment+1:end)],1);
        
        clusterIDvector = zeros(nclusters,1);
        clusterIDvector(Binds) = 1;
        ClusterIDMat = [ClusterIDMat,clusterIDvector];
       
        
    end
    

    
    MetaClusterIDMat(ishuffle,:,:) = ClusterIDMat;

end


SumDot_stim = [];
SumDot_shuffle = [];
for istim = 1:nstims
    
    thisstim = squeeze(MetaClusterIDMat(:,:,istim));
    
    pairind = 1;
    for ishuffle1 = 1:iterations-1
        for ishuffle2 = ishuffle1+1:iterations
            SumDot_shuffle(istim,pairind) = sum(thisstim(ishuffle1,:).*thisstim(ishuffle2,:));
            pairind = pairind+1;
        end
    end
end
for ishuffle = 1:iterations
    
    thisshuffle = squeeze(MetaClusterIDMat(ishuffle,:,:));
    
    pairind = 1;
    for istim1 = 1:nstims-1
        for istim2 = istim1+1:nstims
            SumDot_stim(ishuffle,pairind) = sum(thisshuffle(istim1,:).*thisshuffle(istim2,:));
            pairind = pairind+1;
        end
    end
    
end

estimate_overlap_stims = mean(SumDot_stim,2)
estimate_overlap_shuffles = mean(SumDot_shuffle,2)

[h,p] = kstest2(estimate_overlap_stims,estimate_overlap_shuffles)

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
for ishuffle = 1:iterations
    
    for istim = 1:nstims
%         stimB = abs(squeeze(MetaBetaWeightMatrix(i,:,istim)));
        stimB = (squeeze(MetaBetaWeightMatrix(ishuffle,:,istim)));
        
        ci = 85;
        a = getCDFconf(stimB(find(stimB<0)),ci);
        a = (find(stimB>=a(1)));
        b = getCDFconf(stimB(find(stimB>0)),ci);
        b = (find(stimB<=b(1)));
        signifB_inds = intersect(a,b);

%         signifB_inds = find(stimB);
        
        stimB(signifB_inds) = 0;
        stimB = abs(stimB);
        
%         
%         stimB = abs(stimB(signifB_inds));

        these_locations = allchanxz(include_inds);
%         these_locations = these_locations(signifB_inds);

        shuffleinds = randperm(size(stimB,2));
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
    [x,ishuffle] = sort(stimB);
    
    %     increment = signifCorr_Increment(istim);
    increment = size(corrInc,2);
    
    stimVmHat = squeeze(VmHat(istim,:,:));
    stimVm = Vm(:,istim);
    stimPVm = squeeze(PVm(istim,:,:));
    
    
    Bneg = ishuffle(1:increment);
    %         Bneg = [];
    Bpos = ishuffle(end-increment+1:end);
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
    [x,ishuffle] = sort(stimB);
    
    %     this_increment = round(cell_signifCorr(istim)/100*nclusters)/2;
    %     increment = find(percentTerms == cell_signifCorr(istim));
    increment = signifCorr_Increment(istim);
    if ~isnan(increment)
        Bneg = ishuffle(1:increment);
        %         Bneg = [];
        Bpos = ishuffle(end-increment+1:end);
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

