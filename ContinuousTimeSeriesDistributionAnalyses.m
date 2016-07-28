

%% Extracellular set up cellstruc
cd('/Users/kperks/mnt/cube/Ice/kperks/B1087/MetaMat/')

%from MetaMat folder, where TrialResponsesByStim is a folder containing,
%you guessed it... preprocessed trial responses saved in a separate mat for
%each stimulus

d = dir('TrialResponsesByStim/');
stimmats = [];
for id = 3:size(d,1)
    stimmats{id-2} = d(id).name;
end

load('metatoes.mat','metatoes','xtime_spikes')
fs = metatoes{1}.fs;
stimendtime = metatoes{1}.stims{1}.stim_end_times-metatoes{1}.stims{1}.stim_start_times;
stimdur = stimendtime(1)/fs;

sigon = min(find(xtime_spikes>0));
sigoff = max(find(xtime_spikes<stimdur));

calcstart = sigon + round((sigoff-sigon)/6);

mumean_stim = [];
mumean_prestim = [];
muvar_stim = [];
muvar_prestim = [];
muskew_stim = [];
muskew_prestim = [];
trialvar_stim = [];
trialvar_prestim = [];

for istim = 1:size(stimmats,2)
    load(['TrialResponsesByStim/' stimmats{istim}],'SiteResp')
    
    for isite = 1:size(SiteResp,1)
        if doplot ==1
            figure;hold on
            line(xtime_spikes,squeeze(SiteResp(isite,:,:)),'color','k')
            line(xtime_spikes,mean(squeeze(SiteResp(isite,:,:))),'color','r','LineWidth',4)
            axis tight
            ylabel(['B1083 Site ' num2str(isite) ' Song ' metatoes{1}.stims{1}.name],'Interpreter','none')
        end
    end  
    for isite = 1:size(SiteResp,1)
        thisstim = squeeze(SiteResp(isite,:,calcstart:sigoff));
        thisprestim = squeeze(SiteResp(isite,:,1:calcstart));
        %         for itrial = 1:size(thisstim,1)
        %             stimVm(itrial) = mean(thisstim(itrial,:));
        %             prestimVm(itrial) = mean(thisprestim(itrial,:));
        %         end
        
        mumean_stim(istim,isite,:) = mean(thisstim,2);
        mumean_prestim(istim,isite,:) = mean(thisprestim,2);
        
        muvar_stim(istim,isite,:) = var(thisstim')';
        muvar_prestim(istim,isite,:) = var(thisprestim')';
        
        muskew_stim(istim,isite,:) = skewness(thisstim');
        muskew_prestim(istim,isite,:) = skewness(thisprestim');
        %         clear stimVm prestimVm
        
        trialvar_stim(istim,isite) = mean(var(thisstim));
        trialvar_prestim(istim,isite) = mean(var(thisprestim));
    end
end
%%
colors = colormap(cool(size(trialvar_stim,2)));


hfig = figure;hold on
for icell = 1:size(trialvar_stim,2)
    scatter(trialvar_stim(:,icell),trialvar_prestim(:,icell),200,colors(icell,:))
end
SetAxisUnity(hfig)
% line([0,10],[0,10])
ylabel('prestim  trial(noise) variance PVm')
xlabel('stim period trial(noise) variance PVm')
set(gca,'TickDir','out')


a = mean(mumean_stim,3);
b = mean(mumean_prestim,3);

hfig = figure;hold on
for icell = 1:size(a,2)
    scatter(a(:,icell),b(:,icell),200,colors(icell,:))
end
SetAxisUnity(hfig)
% line([0,10],[0,10])
ylabel('prestim estimate (avg trials) PVm')
xlabel('stim period estimate (avg trials) PVm')
set(gca,'TickDir','out')

aa = mean(muvar_stim,3);
bb = mean(muvar_prestim,3);

hfig = figure;hold on
for icell = 1:size(a,2)
    scatter(aa(:,icell),bb(:,icell),200,colors(icell,:))
end
SetAxisUnity(hfig)
ylabel('prestim tempral variance')
xlabel('stim period temporal variance')
set(gca,'TickDir','out')

% joint distribution of Vm mean and Vm variance for stim and prestim
figure;hold on
scatter(a(:),aa(:),200,'r');
scatter(b(:),bb(:),200,'b');
for icell = 1:size(a,2)
    scatter(nanmean(a(:,icell)),nanmean(aa(:,icell)),200,'filled','r');
    scatter(nanmean(b(:,icell)),nanmean(bb(:,icell)),200,'filled','b');
end
% line([0,35],[0,35])
xlabel('mean Vm')
ylabel('temporal variance')
set(gca,'TickDir','out')
% separated mostly by mean Vm, no information in variance (could do stats on orientation of
% disciminating line for this statement?)

aaa = mean(muskew_stim,3);
bbb = mean(muskew_prestim,3);
hfig = figure;hold on
for icell = 1:size(a,2)
    scatter(aaa(:,icell),bbb(:,icell),200,colors(icell,:))
end
SetAxisUnity(hfig)
ylabel('prestim skewness')
xlabel('stim period skewness')
set(gca,'TickDir','out')

% joint distribution of Vm mean and Vm skewness for stim and prestim

figure;hold on
scatter(aaa(:),a(:),200,'r');
scatter(bbb(:),b(:),200,'b');
for icell = 1:size(a,2)
    scatter(nanmean(aaa(:,icell)),nanmean(a(:,icell)),200,'filled','r');
    scatter(nanmean(bbb(:,icell)),nanmean(b(:,icell)),200,'filled','b');
end
% line([0,35],[0,35])
ylabel('mean Vm')
xlabel('skewness Vm distribution')
set(gca,'TickDir','out')

figure;hold on
for icell = 1:size(a,2)
    scatter(nanmean(aaa(:,icell)),nanmean(a(:,icell)),200,'filled','r');
    scatter(nanmean(bbb(:,icell)),nanmean(b(:,icell)),200,'filled','b');
end
% line([0,35],[0,35])
ylabel('mean Vm')
xlabel('skewness Vm distribution')
% separated mostly by mean Vm, but some information in variance (could do stats on orientation of
% disciminating line for this statement?)

% joint distribution of Vm variance and Vm skewness for stim and prestim

figure;hold on
scatter(aaa(:),aa(:),200,'r');
scatter(bbb(:),bb(:),200,'b');
for icell = 1:size(a,2)
    scatter(nanmean(aaa(:,icell)),nanmean(aa(:,icell)),200,'filled','r');
    scatter(nanmean(bbb(:,icell)),nanmean(bb(:,icell)),200,'filled','b');
end
% line([0,35],[0,35])
ylabel('variance Vm')
xlabel('skewness Vm distribution')
set(gca,'TickDir','out')

%% trained versus novel
%after doing prep on all stim
mumean_stim_N = mumean_stim(1:4,:,:);
muvar_stim_N = muvar_stim(1:4,:,:);
muskew_stim_N = muskew_stim(1:4,:,:);
trialvar_stim_N = trialvar_stim(1:4,:,:);

%after doing prep on all stim
mumean_stim_T = mumean_stim(5:end,:,:);
muvar_stim_T = muvar_stim(5:end,:,:);
muskew_stim_T = muskew_stim(5:end,:,:);
trialvar_stim_T = trialvar_stim(5:end,:,:);

i = mean(mumean_stim_N,3);
j = mean(mumean_stim_T,3);

hfig = figure;hold on
for icell = 1:size(a,2)
    scatter(i(:,icell),j(:,icell),200,colors(icell,:))
end
SetAxisUnity(hfig)
% line([0,10],[0,10])
xlabel('mean PVm novel')
ylabel('mean PVm trained')
set(gca,'TickDir','out')

i = mean(muvar_stim_N,3);
j = mean(muvar_stim_T,3);

hfig = figure;hold on
for icell = 1:size(a,2)
    scatter(i(:,icell),j(:,icell),200,colors(icell,:))
end
SetAxisUnity(hfig)
% line([0,10],[0,10])
xlabel('temporal var PVm novel')
ylabel('temporal PVm trained')
set(gca,'TickDir','out')

i = mean(muskew_stim_N,3);
j = mean(muskew_stim_T,3);

hfig = figure;hold on
for icell = 1:size(a,2)
    scatter(i(:,icell),j(:,icell),200,colors(icell,:))
end
SetAxisUnity(hfig)
% line([0,10],[0,10])
xlabel('skewness PVm novel')
ylabel('skewness PVm trained')
set(gca,'TickDir','out')

i = mean(trialvar_stim_N,3);
j = mean(trialvar_stim_T,3);

hfig = figure;hold on
for icell = 1:size(a,2)
    scatter(i(:,icell),j(:,icell),200,colors(icell,:))
end
SetAxisUnity(hfig)
% line([0,10],[0,10])
xlabel('trial variance PVm novel')
ylabel('trial variance PVm trained')
set(gca,'TickDir','out')
%% estimate of maximal order correlation in synaptic input population


%%
%%%%%%%%%%%%
%INtracellular
%%%%%%%%%%%%%%%%%
BirdList = {'B1075','B1083','B1087','B1235','B1056','st1215'};
a  =[];
aa = [];
aaa = [];
b = [];
bb = [];
bbb = [];
noise_stim = [];
noise_prestim = [];

trainedinds= {[]};

cellind = 1;
for ibird = 1:size(BirdList,2)
    cd(['/Users/kperks/mnt/cube/Ice/kperks/' BirdList{ibird} '/MetaMat/']);
    
    load('IntracellularData.mat')
    % mean of Vm distribution for each cell for pre-stim and stim period
    
    for icell = 1:size(IntracellularData,2)
        mumean_stim = [];
        mumean_prestim = [];
        muvar_stim = [];
        muvar_prestim = [];
        muskew_prestim = [];
        muskew_stim = [];
        muTrialVar_stim = [];
        muTrialVar_prestim = [];
        
        cellstruc = IntracellularData{icell};
        for istim = 1:size(cellstruc.stimnames,2)
            %get Vm distribution per trial
            thisstim = cellstruc.stim_data{istim};
            thisprestim = cellstruc.prestim_data{istim};
            for itrial = 1:size(thisstim,1)
                stimVm(itrial) = mean(thisstim(itrial,:));
                prestimVm(itrial) = mean(thisprestim(itrial,:));
            end
            
            mumean_stim{istim} = (stimVm);
            mumean_prestim{istim} = (prestimVm);
            clear stimVm prestimVm
            
            for itrial = 1:size(thisstim,1)
                stimVm(itrial) = var(thisstim(itrial,:));
                prestimVm(itrial) = var(thisprestim(itrial,:));
            end
            
            muvar_stim{istim} = (stimVm);
            muvar_prestim{istim} = (prestimVm);
            clear stimVm prestimVm
            
            for itrial = 1:size(thisstim,1)
                stimVm(itrial) = skewness(thisstim(itrial,:));
                prestimVm(itrial) = skewness(thisprestim(itrial,:));
            end
            
            muskew_stim{istim} = (stimVm);
            muskew_prestim{istim} = (prestimVm);
            clear stimVm prestimVm
            
            muTrialVar_stim(istim) = mean(var(thisstim));
            muTrialVar_prestim(istim) = mean(var(thisprestim));
        end
    a{cellind} = cellfun(@mean,mumean_stim);
    b{cellind} = cellfun(@mean,mumean_prestim);
    
    aa{cellind} = cellfun(@mean,muvar_stim);
    bb{cellind} = cellfun(@mean,muvar_prestim);
    
    aaa{cellind} = cellfun(@mean,muskew_stim);
    bbb{cellind} = cellfun(@mean,muskew_prestim);
    
    noise_stim{cellind} = muTrialVar_stim;
    noise_prestim{cellind} = muTrialVar_prestim;
    
    cellind = cellind+1;
    end
    
   
    
end
%%
colors = colormap(cool(size(a,2)));


hfig = figure;hold on
for icell = 1:size(a,2)
    scatter(noise_stim{icell}(:),noise_prestim{icell}(:),200,colors(icell,:))
end
SetAxisUnity(hfig)
ylabel('prestim trial(noise) variance Vm')
xlabel('stim trial(noise) variance meanVm')
set(gca,'TickDir','out')


hfig = figure;hold on
for icell = 1:size(a,2)
    scatter(a{icell}(:),b{icell}(:),200,colors(icell,:))
end
SetAxisUnity(hfig)
ylabel('prestim mean Vm')
xlabel('stim period meanVm')
set(gca,'TickDir','out')

hfig = figure;hold on
for icell = 1:size(a,2)
    scatter(aa{icell}(:),bb{icell}(:),200,colors(icell,:))
end
SetAxisUnity(hfig)
ylabel('prestim tempral variance')
xlabel('stim period temporal variance')
set(gca,'TickDir','out')

% joint distribution of Vm mean and Vm variance for stim and prestim

figure;hold on

for icell = 1:size(a,2)
    scatter(a{icell}(:),aa{icell}(:),200,'r');
    scatter(b{icell}(:),bb{icell}(:),200,'b');
    scatter(nanmean(a{icell}(:)),nanmean(aa{icell}(:)),200,'filled','r');
    scatter(nanmean(b{icell}(:)),nanmean(bb{icell}(:)),200,'filled','b');
end
% line([0,35],[0,35])
xlabel('mean Vm')
ylabel('temporal variance')
set(gca,'TickDir','out')

hfig = figure;hold on
for icell = 1:size(a,2)
    scatter(aaa{icell}(:),bbb{icell}(:),200,colors(icell,:))
end
SetAxisUnity(hfig)
ylabel('prestim skewness')
xlabel('stim period skewness')
set(gca,'TickDir','out')

% joint distribution of Vm mean and Vm skewness for stim and prestim

figure;hold on
for icell = 1:size(a,2)
    scatter(aaa{icell}(:),a{icell}(:),200,'r');
    scatter(bbb{icell}(:),b{icell}(:),200,'b');
    scatter(nanmean(aaa{icell}(:)),nanmean(a{icell}(:)),200,'filled','r');
    scatter(nanmean(bbb{icell}(:)),nanmean(b{icell}(:)),200,'filled','b');
end
ylabel('mean Vm')
xlabel('skewness Vm distribution')
set(gca,'TickDir','out')

figure;hold on
for icell = 1:size(a,2)
    scatter(nanmean(aaa{icell}(:)),nanmean(a{icell}(:)),200,'filled','r');
    scatter(nanmean(bbb{icell}(:)),nanmean(b{icell}(:)),200,'filled','b');
end
ylabel('mean Vm')
xlabel('skewness Vm distribution')
% separated mostly by mean Vm, but some information in variance (could do stats on orientation of
% disciminating line for this statement?)

% joint distribution of Vm variance and Vm skewness for stim and prestim

figure;hold on
for icell = 1:size(a,2)
    scatter(aaa{icell}(:),aa{icell}(:),200,'r');
    scatter(bbb{icell}(:),bb{icell}(:),200,'b');
    scatter(nanmean(aaa{icell}(:)),nanmean(aa{icell}(:)),200,'filled','r');
    scatter(nanmean(bbb{icell}(:)),nanmean(bb{icell}(:)),200,'filled','b');
end
ylabel('variance Vm')
xlabel('skewness Vm distribution')
set(gca,'TickDir','out')

figure;hold on
for icell = 1:size(a,2)
    scatter3(aaa{icell}(:),aa{icell}(:),a{icell}(:),200,'r');
    scatter3(bbb{icell}(:),bb{icell}(:),b{icell}(:),200,'b');
    scatter3(nanmean(aaa{icell}(:)),nanmean(aa{icell}(:)),nanmean(a{icell}(:)),200,'filled','r');
    scatter3(nanmean(bbb{icell}(:)),nanmean(bb{icell}(:)),nanmean(b{icell}(:)),200,'filled','b');
end
ylabel('variance Vm')
xlabel('skewness Vm distribution')
zlabel('mean Vm')
set(gca,'TickDir','out')


%% estimate of maximal order correlation in synaptic input population