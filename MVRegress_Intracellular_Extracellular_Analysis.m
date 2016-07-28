

cd('/Users/kperks/mnt/cube/Ice/kperks/B1083/MetaMat/')
currentpath = pwd

path_to_phy = '../klusta/phy040516/';

load('MVRegressPrep.mat','SignalFiltR','include_inds','isNarrow','ExtracellStimNames','IntracellStimNames')

CellID = 1; %which intracellular cell

load('IntracellularData.mat','IntracellularData')
dt = IntracellularData{1}.dt;

load('metatoes.mat','allchanxz','metatoes')
load('behavior.mat')
loc = allchanxz(include_inds);
allinds = [1:size(include_inds,1)];

%%

%%%%%
% also try the regression with pre-prepping wide units with a positive
% kernel and narrow units with a negative kernel
%%%%%
% Pop =[(SignalFiltR);-(SignalFiltR_narrow(:,1,:))];
NarrowInds = isNarrow(include_inds);
Pop = SignalFiltR;
%  Pop = SignalFiltR(find(NarrowInds==0),:,:);
% Pop(find(NarrowInds==1),:,:) = -Pop(find(NarrowInds==1),:,:);
calcstart = size(Pop,3)/6;

% subset_win_Vm = [1:size(Pop,3)-calcstart];
subset_win = [calcstart+1:size(Pop,3)];
% subset_win = [1:calcstart];

nstims = size(SignalFiltR,2);
nclusters = size(SignalFiltR,1);
% nsamps = size(Pop,3);
nsamps = size(subset_win,2);



% if trainonsubset ==0
%     nsamps = size(Pop,3);
% elseif trainonsubset ==1
%     nsamps = size(:,2);
% end


Vm = zeros(nsamps,nstims);
BetaWeightMatrix = zeros(nclusters,nstims);
PVm = zeros(nstims,nsamps,nclusters);
VmHat = zeros(nstims,nsamps,nclusters);
SignalPower = zeros(nsamps,nstims);

pairwiseTrialVar = [];
pairwiseTrialCorr = [];
pairwiseSumR2 = [];

idx = true(2);
idx = ~tril(idx);

for istim = 1:nstims
    istim
    
    %load membrane potential from intracellular cell,
    
    trialavgs = cellfun(@(x) mean(x,1),IntracellularData{CellID}.stim_data,'UniformOutput',0);
    trialavgs = cell2mat(trialavgs');
    globaloffset = min(min(trialavgs));
    trialavgs = trialavgs - globaloffset;
    globalnormfactor = max(max(trialavgs));
    trialavgs = trialavgs ./ globalnormfactor;
    
    Vmdata = mean(IntracellularData{CellID}.stim_data{istim},1);
    Vmdata = Vmdata(subset_win);
    %     Vmdata = Vmdata - min(Vmdata(1,:));
    %     Vmdata = Vmdata ./ max(Vmdata(1,:));
    % for some reason Vmdata is a sample or two longer than Pop data... is
    % a drop in the bucket as far as time, but maybe I can go back and
    % check where it happens and fix it if I end up running the prep again.
    Vmdata = Vmdata - globaloffset;
    Vmdata = Vmdata ./ globalnormfactor;
%     Vmdata = Vmdata(1:size(Pop,3));
%     Vm(:,istim) = fliplr(Vmdata)';
    Vm(:,istim) = Vmdata';
    
    % Pop (SignalFiltR from load) is already "offset to zero" because the
    % traces are built from spiking where no spikes is zero
    % get global normalization factor though for each cluster across stims
    tmpmax = abs((max(max(Pop,[],3),[],2)));
    tmpmin = abs((min(min(Pop,[],3),[],2)));
    cluster_globalnorm = max([tmpmax,tmpmin],[],2);
    
    
    Pop_stim = squeeze(Pop(:,istim,subset_win));
    Pop_stim = Pop_stim';
    
    % for each cluster trial avg, normalize to max for that response on
    % this stimulus (actually this might be messing with asking the
    % question of, across stimuli, are cluster identities consistent.
    %%% I think I should normalize all Vm to stimulus-inclusive max and
    %%% min for that cell (not for the quality metrics calculated across
    %%% trials for a given stim, but for the training data
    PVmdata = [];
    for iterm = 1:size(Pop_stim,2)
        normterm = cluster_globalnorm(iterm) ;
        PVmdata(:,iterm) = Pop_stim(:,iterm)./normterm;
    end
    
    %not sure why some of SignalFiltR would be NaN; maybe it was
    %initialized as that and for some clusters there is sometimes acutally
    %no spikes for a whole stim?
    %for now turn the nan to zero and then investigate later
    PVmdata(isnan(PVmdata)) = 0;
    %     keepterms = find(sum(isnan(PVmdata),1)==0);
    %     PVmdata = PVmdata(:,keepterms);
    
    %     PVm(istim,:,allinds(keepterms)) = PVmdata;
    PVm(istim,:,:) = PVmdata;
    
    [wav,fs] = audioread(['/Users/kperks/mnt/cube/Ice/6sPseudoSongs/' IntracellStimNames{istim} '.wav']);
    y = hilbert(wav);
    env = abs(y);
    dnsampV = round([1:size(env,1)/size(Pop,3):size(env,1)]);
    dnsampWav = env(dnsampV,:);
    SignalPower(:,istim) = dnsampWav(subset_win,1);
end


%%

%%%%%%%%%
%%%% for normal regression
% stiminds = [1:nstims];
%%%% for control/shuffled regression
stiminds = [4:nstims,[1:3]];
%%%%%%%%

for imot = 1%:6-1
%         subset_win = [imot*calcstart:(imot+1)*calcstart];
    % subset_win = [1:size(Vmdata,2)];
    % subset_win = [1:calcstart];%:6*calcstart];
    for istim = 1:nstims
        thisVm = Vm(:,stiminds(istim));
        thisPop = squeeze(PVm(istim,:,:));
        
%         [B] = mvregress(squeeze(PVm(istim,:,:)),Vm(:,stiminds(istim)));
        [B] = regress(thisVm,thisPop);

        BetaWeightMatrix(:,istim) = B;
        
        VmHatdata = [];
        for iterm = 1:size(B,1)
            VmHatdata(:,iterm) = B(iterm)*PVm(istim,:,iterm);
        end
        VmHat(istim,:,:) = VmHatdata;
        
        figure;hold on
        xtime = [1:size(VmHatdata,1)]*dt;
        line(xtime,sum(VmHatdata,2)')
        line(xtime,thisVm,'color','r')
%         legend('VmHat trained only on intro motif','Vm')
        set(gcf,'Position',[298         639        1306         459])
        title(['songID ' num2str(istim)])
    end
end
% figure; hold on
% for istim = 1:nstims
%      xtime = [1:size(VmHat,2)]*dt;
%         line(xtime,sum(squeeze(VmHat(istim,:,:)),2)')
% end
% line(xtime,Vm(:,istim),'color','r')
% set(gcf,'Position',[298         639        1306         459])
%

%%%%%%%%%%%%%%%%%%%%%%%%%%

for istim = 1:nstims
    
    %???? maybe not???? offset and normalize because it makes sense for baseline metrics across trials used for
    % benchmark "quality" to be obtained from comparing across trials the
    % way that PVm and VmHat will be compared to Vm
    thisVm = IntracellularData{CellID}.stim_data{stiminds(istim)};
    
    thisVm = thisVm - repmat(globaloffset,size(thisVm,1),size(thisVm,2));
    thisVm = thisVm ./ repmat(globalnormfactor,size(thisVm,1),size(thisVm,2));
  
    
    %several different metrics for goodness of match to Vm
    pairind = 1;
    for itrial1 = 1:size(thisVm,1)-1
        for itrial2 = itrial1+1:size(thisVm,1)
            
%             pairwiseTrialVar{istim}(pairind) = mean(var(thisVm([itrial1,itrial2],:)));
            
            c = corrcoef(thisVm([itrial1,itrial2],:)');
            pairwiseTrialCorr{istim}(pairind) = c(idx);
            
%             pairwiseSumR2{istim}(pairind) = sum(power(diff(thisVm([itrial1,itrial2],:)),2));
            
            pairind = pairind+1;
        end
    end
    %     figure;   hold on
    %     line(xtime,SignalPower./max(SignalPower))
    %     line(xtime,sum(PVm,2)./max(sum(PVm,2)),'color','m')
    %     line(xtime,Vm./max(Vm),'color','r')
    %     line(xtime,sum(VmHat,2)./max(sum(VmHat,2)),'color','b')
    
    
end
MeanCorrs = cellfun(@(x) mean(x),pairwiseTrialCorr)'
% MeanR2 = round(cellfun(@(x) mean(x),pairwiseSumR2))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

popInc = [];
varInc = [];
corrInc = [];
sumR2Inc = [];
SignalVmCorr = [];
SignalPVmCorr = [];
for istim = 1:nstims
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
%         Bneg = [];
        Bpos = i(end-increment+1:end);
        Binds = [Bneg;Bpos];
        thisB = stimB(Binds,1);
        thisPop = sum(stimVmHat(:,Binds),2);
        %         thisPop = sum(stimVmHat(:,i(end-increment+1:end)),2);
%         thisPop = (thisPop-min(thisPop));
%         thisPop = thisPop./max(thisPop);
        popInc(istim,:,increment) = thisPop;
%         varInc(istim,increment) = mean(var([thisPop';stimVm']));
        c = corrcoef([thisPop';stimVm']');
        corrInc(istim,increment) =c(idx);
%         sumR2Inc(istim,increment) = sum(power(diff([thisPop';stimVm']),2));
    end
    
end
max_corr = max(corrInc')'
% min_R2 = round(min(sumR2Inc'))'
% 
% figure;hold on
% for istim = 1:nstims
%     subplot(nstims/2,2,istim);
%     hold on
% [x, p] = empcdf(corrInc_TrainedIntro(istim,:));          % compute empirical cdf
% stairs(x, p,'color','r');                % plot the cdf
% [x, p] = empcdf(corrInc(istim,:));          % compute empirical cdf
% stairs(x, p);                % plot the cdf
% legend('increment corr Vm VmHat intro','increment corr Vm VmHat whole song')
% set(gca, 'YLim', [0 1]);     % make sure y-range is [0,1]
% end
% 
% figure;hold on
% for istim = 1:nstims
%     subplot(nstims/2,2,istim);
%     hold on
% x = (([1:size(corrInc,2)]*2)/nclusters)*100;
% stairs(x, corrInc_TrainedIntro(istim,:),'color','r');                % plot the cdf
% stairs(x, corrInc(istim,:));          % compute empirical cdf
%            % plot the cdf
% legend('increment corr Vm VmHat intro','increment corr Vm VmHat whole song')
% set(gca, 'YLim', [0 1]);     % make sure y-range is [0,1]
% end
% 
% figure;hold on
% x = (([1:size(corrInc,2)]*2)/nclusters)*100;
% line(x, corrInc(istim,:),'color','k','LineWidth',4);
% line(xlims,[mean(pairwiseTrialCorr{istim}),mean(pairwiseTrialCorr{istim})],'color','r')
% set(gca,'TickDir','out')
% ylabel('correlation coefficient')
% xlabel('percent of terms included in estimate')

%%%%%%%%%%%%%%%%%%%%%
% plot results for Corr metric
percentTerms = (([1:size(sumR2Inc,2)]*2)/nclusters)*100;
nIncrement = [1:size(sumR2Inc,2)];

cell_signifCorr = nan(1,nstims);
% cell_signifSumR2 = nan(1,nstims);
signifCorr_Increment = nan(1,nstims);
for istim = 1:nstims
    
    if ~isempty(find(corrInc(istim,:)<mean(pairwiseTrialCorr{istim})))
        cell_signifCorr(istim) = percentTerms(min(find(corrInc(istim,:)>mean(pairwiseTrialCorr{istim}))));
        signifCorr_Increment(istim) = nIncrement(min(find(corrInc(istim,:)>mean(pairwiseTrialCorr{istim}))));
    else 
        cell_signifCorr(istim) = 2/nclusters*100;
        signifCorr_Increment(istim) = 1;
    end
    
    % it doesn't seem to make sense to me to use R2 unless I am normalizing
    % every trace that i use in all of these analyses (which also does not
    % make sense)
    % AND when i do normalize traces and use R2 the result seems consistent
    % with using corr as the assessment metric 
%     if ~isempty(find(sumR2Inc(istim,:)<mean(pairwiseSumR2{istim})))
%         cell_signifSumR2(istim) = percentTerms(min(find(sumR2Inc(istim,:)<mean(pairwiseSumR2{istim}))));
%     else
%         cell_signifSumR2(istim) = 2/nclusters*100;
%     end
end
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

%% Wide / Narrow mean population level activity of only predictive untis
for istim = 1:nstims
    stimB = BetaWeightMatrix(:,istim);
    [x,i] = sort(stimB);
    
    increment = signifCorr_Increment(istim);
%     increment = size(corrInc,2);
    
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
title(['stimulus ' num2str(istim) ' B1083 n=' num2str(increment*2) ' units in population was necessary'])

echo_corrBnegBpos = corrcoef(thisPopHat_Bneg,thisPopHat_Bpos)
echo_corrBnegVm = corrcoef(thisPopHat_Bneg,stimVm)
echo_corrBposVm = corrcoef(thisPopHat_Bpos,stimVm)
echo_corrSumPopVm = corrcoef(sum([thisPopHat_Bpos,thisPopHat_Bneg],2),stimVm)
%% Narrow is negative Beta weights?
NarrowBool = find(NarrowInds==1);
WideBool = find(NarrowInds==0);
nNarrow = size(NarrowBool,2);
nWide = size(WideBool,2);

HitN = zeros(1,nstims);
HitP = zeros(1,nstims);
FaP = zeros(1,nstims);
FaN = zeros(1,nstims);
for istim = 1:nstims
    stimB = BetaWeightMatrix(:,istim);
    [Bsorted,i_order] = sort(stimB);
    
    isBneg = i_order(find(Bsorted<0));
    isBpos = i_order(find(Bsorted>0));

    nBneg = size(isBneg,1);
    nBpos = size(isBpos,1);
    
    for i = 1:nNarrow
        if ~isempty(find(isBneg == NarrowBool(i)))
            HitN(istim) = HitN(istim)+1;
        end
        if ~isempty(find(isBpos == NarrowBool(i)))
            FaP(istim) = FaP(istim)+1;
        end
        
    end
    for i = 1:nWide
        if ~isempty(find(isBpos == WideBool(i)))
            HitP(istim) = HitP(istim)+1;
        end
        if ~isempty(find(isBneg == WideBool(i)))
            FaN(istim) = FaN(istim)+1;
        end
        
    end
    
end

AccNarrow = HitN ./ (HitN + FaP);
AccWide = HitP ./ (HitP + FaN);
figure;hold on
scatter(AccNarrow,AccWide,200,'fill','k')
ylabel('accuracy wide categorized w/ pos B')
xlabel('accuaracy narrow categorized w/ neg B')
legend('n = 8 songs')

figure;hold on
scatter(repmat(1,1,nstims),HitN,200,'fill','r')
scatter(repmat(2,1,nstims),FaN,200,'r')
scatter(repmat(3,1,nstims),HitP,200,'fill','b')
scatter(repmat(4,1,nstims),FaP,200,'b')
set(gca,'XLim',[0,5],'XTick',[1:4],'XTickLabel',{'Hit Neg','FA Neg','Hit Pos','FA Pos'})

%% plot results for R2
location_mat = zeros(size(unique(loc),1),nstims);
loc_array = unique(loc);
for istim = 1:nstims
    stimB = BetaWeightMatrix(:,istim);
    [x,i] = sort(stimB);
    
    this_increment = round(cell_signifSumR2(istim)/100*nclusters)/2;
    
    if ~isnan(this_increment)
        these_inds = i([i(1:this_increment),i(end-this_increment+1:end)],1);
        these_locations = loc(these_inds);
        for ib = 1:size(these_inds,1)
            location_mat(find(loc_array == these_locations(ib)),istim) = 1;
        end
    end
end
figure;imagesc(location_mat)
set(gca,'YTick',[1:4:size(loc_array,1)],'YTickLabel',loc_array([1:4:size(loc_array,1)]))
set(gca,'XTick',[1:nstims],'XTickLabel',metatoes_songlist)
ylabel('location along penetration')
xlabel('PseudoSong label')
title('significance by R2')
colorbar

figure;hold on
percentTerms = (([1:size(sumR2Inc,2)]*2)/nclusters)*100;
for istim = 1:nstims
    subplot(nstims/2,2,istim)
    ylabel(['sumR2 stim ' num2str(istim)])
    line(percentTerms,sumR2Inc(istim,:))
    xlims = get(gca,'XLim');
    line(xlims,[mean(pairwiseSumR2{istim}),mean(pairwiseSumR2{istim})],'color','r')
    %     set(gca,'YLim',[0,1])
    
end

istim = 1;
figure; hold on
% [1:size(popInc,2)]
line(subset_win*dt,popInc(istim,:,round(size(varInc,2)*0.05)),'color','m')
line(subset_win*dt,popInc(istim,:,round(size(varInc,2)*0.25)),'color','r')
line(subset_win*dt,popInc(istim,:,round(size(varInc,2)*1)),'color','b')
line(subset_win*dt,Vm(:,istim),'color','k')
set(gca,'YLim',[-1,2])
legend('5% VmHat','25% VmHat','100% VmHat','Vm')
title(['stim ' num2str(istim)])
axis tight
set(gca,'TickDir','out','XTick',[1:6],'YTick',[0,1])

% hfig = figure;hold on
% for icell = 1%:size(signifVar,1)
% scatter(cell_signifCorr(1:4),cell_signifCorr(5:8),200,'fill','r')
% scatter(cell_signifSumR2(1:4),cell_signifSumR2(5:8),200,'fill','b')
% end
% SetAxisUnity(hfig)
% legend('Corr','SumR2')
% ylabel('stims 5:8')
% xlabel('stims 1:4')



%%
Cellind = 3;
percentTerms = (([1:size(sumR2Inc,2)]*2)/nclusters)*100;
signifCorr = [];
signifSumR2 = [];
for istim = 1:nstims
    signifCorr(Cellind,istim) = percentTerms(max(find(corrInc(istim,:)<mean(pairwiseTrialCorr{istim}))));
    signifSumR2(Cellind,istim) = percentTerms(min(find(sumR2Inc(istim,:)<mean(pairwiseSumR2{istim}))));
end

hfig = figure;hold on
for icell = Cellind%:size(signifVar,1)
    scatter(signifCorr(icell,1:4),signifCorr(icell,5:8),200,'fill','r')
    scatter(signifSumR2(icell,1:4),signifSumR2(icell,5:8),200,'fill','b')
end
SetAxisUnity(hfig)
legend('corr','SumR2')
ylabel('stims 5:8')
xlabel('stims 1:4')
%%
%  figure;hold on
% for istim = 1:nstims
%     subplot(nstims/2,2,istim)
%     ylabel(['var stim ' num2str(istim)])
%     line((([1:size(varInc,2)]*2)/nclusters)*100,varInc(istim,:))
%     xlims = get(gca,'XLim');
%     line(xlims,[mean(pairwiseTrialVar{istim}),mean(pairwiseTrialVar{istim})],'color','r')
%     set(gca,'YLim',[0,0.04])
% end

cell_signifCorr = [];
figure;hold on
for istim = 1:nstims
    subplot(nstims/2,2,istim)
    ylabel(['corr stim ' num2str(istim)])
    line((([1:size(corrInc,2)]*2)/nclusters)*100,corrInc(istim,:))
    xlims = get(gca,'XLim');
    line(xlims,[mean(pairwiseTrialCorr{istim}),mean(pairwiseTrialCorr{istim})],'color','r')
    set(gca,'YLim',[0,1])
    
    if ~isempty(find(corrInc(istim,:)<mean(pairwiseTrialCorr{istim})))
        cell_signifCorr(istim) = percentTerms(max(find(corrInc(istim,:)<mean(pairwiseTrialCorr{istim}))));
    end
end

cell_signifSumR2 = [];
figure;hold on
percentTerms = (([1:size(sumR2Inc,2)]*2)/nclusters)*100;
for istim = 1:nstims
    subplot(nstims/2,2,istim)
    ylabel(['sumR2 stim ' num2str(istim)])
    line(percentTerms,sumR2Inc(istim,:))
    xlims = get(gca,'XLim');
    line(xlims,[mean(pairwiseSumR2{istim}),mean(pairwiseSumR2{istim})],'color','r')
    %     set(gca,'YLim',[0,1])
    
    if ~isempty(find(sumR2Inc(istim,:)<mean(pairwiseSumR2{istim})))
        cell_signifSumR2(istim) = percentTerms(min(find(sumR2Inc(istim,:)<mean(pairwiseSumR2{istim}))));
    end
end


hfig = figure;hold on
for icell = Cellind%:size(signifVar,1)
    scatter(cell_signifCorr(1:4),cell_signifCorr(5:8),200,'fill','r')
    scatter(cell_signifSumR2(1:4),cell_signifSumR2(5:8),200,'fill','b')
end
SetAxisUnity(hfig)
legend('Corr','SumR2')
ylabel('stims 5:8')
xlabel('stims 1:4')

istim = 1;
figure; hold on
line([1:nsamps]*dt,popInc(istim,:,round(size(varInc,2)*0.05)),'color','m')
line([1:nsamps]*dt,popInc(istim,:,round(size(varInc,2)*0.30)),'color','r')
line([1:nsamps]*dt,popInc(istim,:,round(size(varInc,2)*1)),'color','b')
line([1:nsamps]*dt,Vm(:,istim),'color','k')
legend('5% VmHat','30% VmHat','100% VmHat','Vm')
title(['stim ' num2str(istim)])

%
% % keep increasing number of beta included until VmHat is within reasonable variance (constrained by empirical trial-trial corr or error)
%     VmHat_signif = [];
%     for iterm = 1:size(signifinds,1)
%         VmHat_signif{istim}(:,iterm) = B(signifinds(iterm))*PVm(:,signifinds(iterm));
%     end
%%

%%%%%
% also try the regression with pre-prepping wide units with a positive
% kernel and narrow units with a negative kernel
%%%%%
% Pop =[(SignalFiltR);-(SignalFiltR_narrow(:,1,:))];
NarrowInds = isNarrow(include_inds);
Pop = SignalFiltR;
%  Pop = SignalFiltR(find(NarrowInds==0),:,:);
Pop(find(NarrowInds==1),:,:) = -Pop(find(NarrowInds==1),:,:);

Beta = [];
Beta_locs = [];
nstims = size(SignalFiltR,2);
for istim = 1:nstims
    istim
    Vm = mean(IntracellularData{CellID}.stim_data{istim},1);
    Vm = Vm - min(Vm);
    %     Vm = medfilt1(Vm,medsmoothing,[],2);
    Vm = Vm(1:size(Pop,3));
    Vm = Vm';
    
    Pop_stim = squeeze(Pop(:,istim,:));
    Pop_stim = Pop_stim';
    %     keepterms = find(sum(isnan(Pop_stim),1)==0);
    PVm = [];
    for iterm = 1:size(Pop_stim,2)
        normterm = max([abs(max(Pop_stim(:,iterm))),abs(min(Pop_stim(:,iterm)))]);
        
        PVm(:,iterm) = Pop_stim(:,iterm)./normterm;
        
    end
    keepterms = find(sum(isnan(PVm),1)==0);
    PVm = PVm(:,keepterms);
    
    [B] = mvregress(PVm,Vm);
    Beta{istim} = B;
    Beta_locs{istim} = loc(keepterms);
    
    if doplot ==1
        figure;hold on
        set(gcf,'Position',[342          80        1538        1018])
        subplot(3,2,1)
        hold on
        xtime=[1:size(Pop_stim,1)]*IntracellularData{1}.dt;
        line(xtime,mean(Pop_stim(:,find(NarrowInds==1)),2),'color','r')
        line(xtime,mean(Pop_stim(:,find(NarrowInds==0)),2),'color','b')
        axis tight
        legend('pop mean narrow events','pop mean wide events')
        
        subplot(3,2,3)
        line(xtime,Vm)
        ylabel('trial avg Intracell Vm')
        
        [wav,fs] = audioread(['/Users/kperks/mnt/cube/Ice/6sPseudoSongs/' IntracellStimNames{istim} '.wav']);
        y = hilbert(wav);
        env = abs(y);
        dnsampV = round([1:size(env,1)/size(VmHat,1):size(env,1)]);
        SignalPower = env(dnsampV,:);
        subplot(3,2,5)
        line([1:length(wav)]/fs,wav)
        axis tight
        ylabel(['Pseudo Song ' IntracellStimNames{istim}]);
        
        B = Beta{istim};
        y = prctile(B, [5 50 95]);
        % signifinds = union(find(B>=y(3)),find(B<=y(1)));
        signifinds = find(B>=y(3));
        
        VmHat = [];
        for iterm = 1:size(B,1)
            VmHat(:,iterm) = B(iterm)*PVm(:,iterm);
        end
        VmHat_signif = [];
        for iterm = 1:size(signifinds,1)
            VmHat_signif(:,iterm) = B(signifinds(iterm))*PVm(:,signifinds(iterm));
        end
        
        %     figure;   hold on
        %     line(xtime,SignalPower./max(SignalPower))
        %     line(xtime,sum(PVm,2)./max(sum(PVm,2)),'color','m')
        %     line(xtime,Vm./max(Vm),'color','r')
        %     line(xtime,sum(VmHat,2)./max(sum(VmHat,2)),'color','b')
        
        subplot(3,2,2);
        hold on
        line(xtime,sum(PVm,2)./max(sum(PVm,2)),'color','k')
        line(xtime,sum(VmHat,2)./max(sum(VmHat,2)),'color','m')
        axis tight
        ylabel('normalized amplitude')
        legend('sum pop','sum weighted pop')
        
        subplot(3,2,4)
        line(xtime,sum(VmHat_signif,2)./max(sum(VmHat_signif,2)),'color','k')
        ylabel('sum VmHat +5% signif')
        axis tight
        
        subplot(3,2,6)
        line(xtime,Vm./max(Vm),'color','k')
        line(xtime,sum(VmHat,2)./max(sum(VmHat,2)),'color','m')
        ylabel('normalized amplitude')
        axis tight
        legend('Vm','sum weighted pop')
    end
    
end

%across a sliding location window, get mean B
loc_breadth = 50;
loc_centers = [min(unique(loc))+loc_breadth:50:max(unique(loc))-loc_breadth];

B_globalmean = [];
Bhat = [];
hfigT = figure;hold on
xlabel('Beta values within 100microns of location')
ylabel('electrode z primary location')
title('trained songs')
hfigN = figure;hold on
xlabel('Beta values within 100microns of each location')
ylabel('electrode z primary location')
title('novel songs')
hfigLocN = figure;hold on
xlabel('n clusters within 100microns of location')
ylabel('elextrode z primary location')
for iB = 1:size(Beta,2)
    
    
    B = Beta{iB};
    
    Bloc = Beta_locs{iB};
    
    for icenter = 1:size(loc_centers,2)
        if ~isempty(find(trainedSong_ind==iB))
            figure(hfigT)
            colormap = 'm';
        elseif ~isempty(find(novelSong_ind==iB))
            figure(hfigN)
            colormap = 'k';
        end
        centerinds = intersect(find(Bloc>loc_centers(icenter)-loc_breadth),find(Bloc<loc_centers(icenter)+loc_breadth));
        theseinds = centerinds;
        %         theseinds = intersect(centerinds, posBinds);
        Bhat(icenter,iB) = mean(B(theseinds));
        scatter(B(theseinds),repmat(icenter,size(theseinds,1),1),100,colormap)
        
        figure(hfigLocN)
        scatter(size(theseinds,1),icenter,200,'k','fill')
    end
    posBinds = find(B>0);
    B_globalmean(iB) = mean(B(posBinds));
    
    
end
figure(hfigT)
ytickinds = get(gca,'YTick');
set(gca,'YTick',ytickinds(2:end-1))
set(gca,'YTickLabel',loc_centers(ytickinds(2:end-1)))
figure(hfigN)
ytickinds = get(gca,'YTick');
set(gca,'YTick',ytickinds(2:end-1))
set(gca,'YTickLabel',loc_centers(ytickinds(2:end-1)))
figure(hfigLocN)
ytickinds = get(gca,'YTick');
set(gca,'YTick',ytickinds(2:end-1))
set(gca,'YTickLabel',loc_centers(ytickinds(2:end-1)))

%across a sliding location window, get mean positive weight B
B_globalmean = [];
Bhat = [];
for iB = 1:size(Beta,2)
    B = Beta{iB};
    posBinds = find(B>0);
    Bloc = Beta_locs{iB};
    
    for icenter = 1:size(loc_centers,2)
        centerinds = intersect(find(Bloc>loc_centers(icenter)-loc_breadth),find(Bloc<loc_centers(icenter)+loc_breadth));
        theseinds = intersect(centerinds, posBinds);
        Bhat(icenter,iB) = mean(B(theseinds));
    end
    B_globalmean(iB) = mean(B(posBinds));
end
figure;hold on
line(loc_centers,mean(Bhat,2),'LineWidth',4)
line([min(loc_centers),max(loc_centers)],[mean(B_globalmean),mean(B_globalmean)],'color',[0.5 0.5 0.5])
ylabel('beta estimate within each location center')
title('estimate (avg across stim) Beta (positive only) within 100microns of center location')


B_globalmean = [];
Bhat = [];
for iB = 1:size(Beta,2)
    B = Beta{iB};
    posBinds = find(B<0);
    Bloc = Beta_locs{iB};
    
    for icenter = 1:size(loc_centers,2)
        centerinds = intersect(find(Bloc>loc_centers(icenter)-loc_breadth),find(Bloc<loc_centers(icenter)+loc_breadth));
        theseinds = intersect(centerinds, posBinds);
        Bhat(icenter,iB) = mean(B(theseinds));
    end
    B_globalmean(iB) = mean(B(posBinds));
end
figure;hold on
line(loc_centers,mean(Bhat,2),'LineWidth',4)
line([min(loc_centers),max(loc_centers)],[mean(B_globalmean),mean(B_globalmean)],'color',[0.5 0.5 0.5])
ylabel('beta estimate within each location center')
title('estimate (avg across stim) Beta (negative only) within 100microns of center location')

%%
titstr = '5% confidence intervals on B for Vhat';

edges = [min(unique(allchanxz)):60:max(unique(allchanxz))];
n_all = zeros(size(edges,2),size(Beta,2));
n_signif = zeros(size(edges,2),size(Beta,2));
for iB = 1:size(Beta,2)
    B = Beta{iB};
    y = prctile(B, [5 50 95]);
    % signifinds = union(find(B>=y(3)),find(B<=y(1)));
    signifinds = find(B>=y(3));
    
    Bloc = Beta_locs{iB};
    
    [n,bins] = histc(Bloc,edges);
    n_all(:,iB) = n./size(Bloc,1);
    
    [n,bins] = histc(Bloc(signifinds),edges);
    n_signif(:,iB) = n./size(signifinds,1);
    
    %     smear = 5;
    %     t = medfilt1(n_signif,smear,[],1);
    %     tt = medfilt1(n_all,smear,[],1);
    %     figure;hold on;
    %     line(edges,tt./max(tt),'color','k')
    %     line(edges,t./max(t),'color','r')
end
figure;hold on
line(edges, mean(n_signif(:,[1:4,9:12]),2),'color','b')
line(edges, mean(n_signif(:,[5:8,13:16]),2),'color','g')
line(edges, n_all,'color','k')
axis tight
legend('mean location distribution for significant beta novel stim','mean location distribution for significant beta trained stim','distribution locations all')

figure;hold on
line(edges, mean(n_signif,2),'color','r')
line(edges, mean(n_all,2),'color','k')
axis tight
legend('mean location distribution for significant beta','distribution locations all')

figure;hold on
line(edges, mean(n_signif,2)./mean(n_all,2),'color','k')
axis tight
title('mean location distribution significant beta ./ all')

%%%%plot for one stim
% istim = last stim in for loop above
figure;scatter([1:size(B,1)],sort(B))
line([1,size(B,1)],[y(1),y(1)])
line([1,size(B,1)],[y(3),y(3)])
axis tight
title(titstr)

VmHat = [];
for iterm = 1:size(B,1)
    VmHat(:,iterm) = B(iterm)*PVm(:,iterm);
end

VmHat_signif = [];
for iterm = 1:size(signifinds,1)
    VmHat_signif(:,iterm) = B(signifinds(iterm))*PVm(:,signifinds(iterm));
end

figure;
hold on
plot(sum(PVm,2)./max(sum(PVm,2)),'color','m')
plot((Vm)./max(Vm),'color','k')
plot(sum(VmHat,2)./max(sum(VmHat,2)),'color','b')
plot(sum(VmHat_signif,2)./max(sum(VmHat_signif,2)),'color','r')
legend('sum PVm','Vm','sum VmHat','sum VmHat signif')
axis tight
title(titstr)
%% correlation between sum population activity and intracellular Vm... novel versus trained songs, does intracellular input reflect more "common field" in one case versus other?
CellID = 2;
load('behavior.mat')
NarrowInds = isNarrow(include_inds);
Pop = SignalFiltR;
%  Pop = SignalFiltR(find(NarrowInds==0),:,:);
Pop(find(NarrowInds==1),:,:) = -Pop(find(NarrowInds==1),:,:);

CorrVmPVm = [];
LagVmPm = [];
Intro_CorrVmPVm = [];
Intro_LagVmPm = [];

calcstart = round(size(Pop,3)/6);

nstims = size(SignalFiltR,2);
hfig = figure;hold on
for istim = 1:nstims
    istim
    Vm = mean(IntracellularData{CellID}.stim_data{istim},1);
    Vm = Vm - min(Vm);
    %     Vm = medfilt1(Vm,medsmoothing,[],2);
    Vm = Vm(1:size(Pop,3));
    Vm = Vm ./ max(Vm);
    Vm_resp = Vm(calcstart:end);
    Vm_intro = Vm(1:calcstart);
    %     Vm = Vm';
    
    Pop_stim = squeeze(Pop(:,istim,:));
    
    %     figure;
    %     hold on;
    %     norm = max(-mean(Pop_stim(find(NarrowInds==1),:)));
    %     offset = min(-mean(Pop_stim(find(NarrowInds==1),:)));
    %     line([1:size(Pop_stim,2)]*IntracellularData{1}.dt,(-mean(Pop_stim(find(NarrowInds==1),:))/norm)-offset,'color','r')
    %     norm = max(mean(Pop_stim(find(NarrowInds==0),:)));
    %     offset = min(mean(Pop_stim(find(NarrowInds==0),:)));
    %     line([1:size(Pop_stim,2)]*IntracellularData{1}.dt,(mean(Pop_stim(find(NarrowInds==0),:))/norm)-offset,'color','b')
    
    PVm = sum(Pop_stim,1);
    PVm = PVm ./ max(PVm);
    PVm_resp = PVm(calcstart:end);
    PVm_intro = PVm(1:calcstart);
    
    [r,lags] = xcorr(PVm_resp,Vm_resp,10000,'coeff');
    CorrVmPm(istim) = unique(max(r));
    LagVmPm(istim) = lags(min(find(r==max(r))));
    
    [r,lags] = xcorr(PVm_intro,Vm_intro,10000,'coeff');
    Intro_CorrVmPm(istim) = unique(max(r));
    Intro_LagVmPm(istim) = lags(min(find(r==max(r))));
    
    subplot(4,2,istim);
    hold on;
    line([1:size(PVm_resp,2)]*IntracellularData{1}.dt,PVm_resp','color','m')
    line([1:size(Vm_resp,2)]*IntracellularData{1}.dt,Vm_resp','color',[0.7,0.5,0])
    axis tight
    ylabel(IntracellularData{CellID}.stimnames{istim})
    %     Pop_stim = Pop_stim';
    %     keepterms = find(sum(isnan(Pop_stim),1)==0);
    %     PVm = [];
    %     for iterm = 1:size(Pop_stim,2)
    %         normterm = max([abs(max(Pop_stim(:,iterm))),abs(min(Pop_stim(:,iterm)))]);
    %
    %         PVm(:,iterm) = Pop_stim(:,iterm)./normterm;
    %
    %     end
    %     keepterms = find(sum(isnan(PVm),1)==0);
    %     PVm = PVm(:,keepterms);
end

axis tight
legend('PVm','Vm')
ylabel('normalized to max')
xlabel('seconds')

figure;hold on
scatter(repmat(2,size(trainedSong_ind,2),1),CorrVmPm(trainedSong_ind),200,'r','fill')
scatter(repmat(1,size(novelSong_ind,2),1),CorrVmPm(novelSong_ind),200,'b','fill')
set(gca,'XTick',[1,2],'XLim',[0,3],'XTickLabel',{'Novel','Trained'},'YLim',[0.5,1])
ylabel('Max Correlation')
title('post-intro response')

figure;hold on
scatter(repmat(2,size(trainedSong_ind,2),1),LagVmPm(trainedSong_ind),200,'r','fill')
scatter(repmat(1,size(novelSong_ind,2),1),LagVmPm(novelSong_ind),200,'b','fill')
set(gca,'XTick',[1,2],'XLim',[0,3],'XTickLabel',{'Novel','Trained'})
ylabel('Lag at Max Correlation (PVm vs Vm)')
title('post-intro response')

figure;hold on
scatter(repmat(2,size(trainedSong_ind,2),1),Intro_CorrVmPm(trainedSong_ind),200,'r','fill')
scatter(repmat(1,size(novelSong_ind,2),1),Intro_CorrVmPm(novelSong_ind),200,'b','fill')
set(gca,'XTick',[1,2],'XLim',[0,3],'XTickLabel',{'Novel','Trained'},'YLim',[0.5,1])
ylabel('Max Correlation')
title('intro response')

figure;hold on
scatter(repmat(2,size(trainedSong_ind,2),1),Intro_LagVmPm(trainedSong_ind),200,'r','fill')
scatter(repmat(1,size(novelSong_ind,2),1),Intro_LagVmPm(novelSong_ind),200,'b','fill')
set(gca,'XTick',[1,2],'XLim',[0,3],'XTickLabel',{'Novel','Trained'})
ylabel('Lag at Max Correlation (PVm vs Vm)')
title('intro response')

%% correlation between wide and narrow sum population activity
fs = 1/IntracellularData{1}.dt;
tau = 0.01;
load('metatoes.mat','xtime_spikes')

NarrowInds = isNarrow(include_inds);
Pop = SignalFiltR;
calcstart = round(size(Pop,3)/6);
%  Pop = SignalFiltR(find(NarrowInds==0),:,:);
Pop(find(NarrowInds==1),:,:) = -Pop(find(NarrowInds==1),:,:);

W_mag = [];
N_mag = [];
W_mag_intro = [];
N_mag_intro = [];
WN_corr = [];
WN_lag = [];

nstims = size(SignalFiltR,2);
for istim = 1:nstims
    
    Pop_stim = squeeze(Pop(:,istim,:));
    Pop_sum_W = sum(Pop_stim(find(NarrowInds==0),:),1);
    Pop_sum_N = sum(Pop_stim(find(NarrowInds==1),:),1);
    
    [r,lags] = xcorr(Pop_sum_W,Pop_sum_N,10000,'coeff');
    WN_corr(istim) = min(r);
    WN_lag(istim) = lags(find(r == min(r)))/fs;
    
    W_mag(:,istim) = (sum(Pop_stim(find(NarrowInds==0),calcstart:end),2)/(tau*fs))/((size(Pop_stim,2)-calcstart)/fs);
    N_mag(:,istim) = (-sum(Pop_stim(find(NarrowInds==1),calcstart:end),2)/(tau*fs))/((size(Pop_stim,2)-calcstart)/fs);
    
    W_mag_intro(:,istim) = (sum(Pop_stim(find(NarrowInds==0),1:calcstart),2)/(tau*fs))/((calcstart)/fs);
    N_mag_intro(:,istim) = (-sum(Pop_stim(find(NarrowInds==1),1:calcstart),2)/(tau*fs))/((calcstart)/fs);
    
    
end

edges = [1:1:round(max([max(max(W_mag)),max(max(N_mag))]))];
nW = histc(W_mag(:),edges);
nW = nW/size(W_mag(:),1);
nN = histc(N_mag(:),edges);
nN = nN/size(N_mag(:),1);
figure;hold on
stairs(edges,nW,'color','b','LineWidth',4)
stairs(edges,nN,'color','r','LineWidth',4)

edges = [0:1:round(max([max(max(W_mag_intro)),max(max(N_mag_intro))]))];
nW_intro = histc(W_mag_intro(:),edges);
nW_intro = nW_intro/size(W_mag_intro(:),1);
nN_intro = histc(N_mag_intro(:),edges);
nN_intro = nN_intro/size(N_mag_intro(:),1);
figure;hold on
stairs(edges,nW_intro,'color','b','LineWidth',4)
stairs(edges,nN_intro,'color','r','LineWidth',4)


figure;scatter(W_mag(:),W_mag_intro(:))
figure;scatter(N_mag(:),N_mag_intro(:))

%% pairwise signal correlation and geometric mean spike-count(ish) for wide and narrow and mixed pairs

NarrowInds = isNarrow(include_inds);
Pop = SignalFiltR;
calcstart = round(size(Pop,3)/6);
%  Pop = SignalFiltR(find(NarrowInds==0),:,:);
% Pop(find(NarrowInds==1),:,:) = -Pop(find(NarrowInds==1),:,:);

nclusters = size(SignalFiltR,1);
nstims = size(SignalFiltR,2);

g = nan(nclusters-1,nclusters,nstims);
c = nan(nclusters-1,nclusters,nstims);
pID = nan(nclusters-1,nclusters);

idx = true(2);
idx = ~tril(idx);

for icluster1 = 1:nclusters-1
    
    for icluster2 = icluster1+1:nclusters
        
        if NarrowInds(icluster1)==0 && NarrowInds(icluster2)==0
            pID(icluster1,icluster2) = 1;
        elseif NarrowInds(icluster1)==1 && NarrowInds(icluster2)==1
            pID(icluster1,icluster2) = 2;
        elseif NarrowInds(icluster1)==0 && NarrowInds(icluster2)==1
            pID(icluster1,icluster2) = 3;
        elseif NarrowInds(icluster1)==1 && NarrowInds(icluster2)==0
            pID(icluster1,icluster2) = 3;
        end
        
        for istim = 1:nstims
            
            R1 = squeeze(Pop(icluster1,istim,calcstart:end));
            R2 = squeeze(Pop(icluster2,istim,calcstart:end));
            c_tmp = corrcoef(R1,R2);
            c(icluster1,icluster2,istim) = c_tmp(idx);
            
            g(icluster1,icluster2,istim) = (sum(R1)/(tau*fs))* (sum(R2)/(tau*fs));
            
            
        end
    end
end

c = c(:,2:end,:);
g = g(:,2:end,:);
pID = pID(:,2:end,:);

idx = true(size(c,1));
% idd = triu(idx).*tril(idx);
% pID(find(idd)) = 0;
idx = ~tril(idx);

c_hat = mean(c,3);
c_hat = c_hat(idx);
g_hat = mean(g,3);
WWcorrs = c_hat(find(pID==1));
NNcorrs = c_hat(find(pID==2));
WNcorrs = c_hat(find(pID==3));

WWgeomean = g_hat(find(pID==1));
NNgeomean = g_hat(find(pID==2));
WNgeomean = g_hat(find(pID==3));

edges = [-1:0.1:1];
figure;hold on
n = histc(WWcorrs,edges);
stairs(edges,n/size(WWcorrs,1),'color','b','LineWidth',4);
n = histc(NNcorrs,edges);
stairs(edges,n/size(NNcorrs,1),'color','r','LineWidth',4);
n = histc(WNcorrs,edges);
stairs(edges,n/size(WNcorrs,1),'color','k','LineWidth',4);
legend('Wide:Wide correlation','Narrow:Narrow corrlation','Wide:Narrow correlation')


edges = [0:1:200];
figure;hold on
n = histc(sqrt(WWgeomean),edges);
stairs(edges,n/size(WWgeomean,1),'color','b','LineWidth',4);
n = histc(sqrt(NNgeomean),edges);
stairs(edges,n/size(NNgeomean,1),'color','r','LineWidth',4);
n = histc(sqrt(WNgeomean),edges);
stairs(edges,n/size(WNgeomean,1),'color','k','LineWidth',4);
legend('Wide:Wide geomean','Narrow:Narrow geomean','Wide:Narrow geomean')

figure;
scatter(sqrt(WWgeomean),WWcorrs)
xlabel('WW geomean')
ylabel('WW pairwise corrs')

figure;
scatter(sqrt(NNgeomean),NNcorrs)
xlabel('NN geomean')
ylabel('NN pairwise corrs')

figure;
scatter(sqrt(WNgeomean),WNcorrs)
xlabel('WN geomean')
ylabel('WN pairwise corrs')
