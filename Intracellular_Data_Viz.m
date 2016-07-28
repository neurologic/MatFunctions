set(0,'DefaultAxesFontSize',20)
set(0,'DefaultTextFontSize',20)
set(0,'DefaultTextFontName','Helvetica')
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesColorOrder',[0 0 0;1 0 1;0 1 1;0.6 0.6 0])


%% see if mean VM correlates with trial #
trialVm = nan(1,size(expt.wc.data,1));
trialBins = expt.sweeps.trial; %[1:size(expt.wc.data,1)];
for itrial = 1:size(trialBins,1)
   thistrial = filtesweeps(expt,0,'trial',trialBins(itrial)) ;
   trialVm(itrial) = mean(thistrial.wc.data);
end
figure;
scatter(trialBins,trialVm);
xlabel('trial number')
ylabel('mean Vm')

title(num2str(unique(filtexpt.sweeps.Vm)))

set(gcf,'Position',[585   894   809   212])

%% filtexpt by trial if needed
filtexpt = filtesweeps(expt,0,'Vm',[0]);
filtexpt = filtesweeps(expt,0,'trial',[6:72]);

%% get spike thresh
highpassdata=HighpassGeneral(filtexpt.wc.data*1000,1/expt.wc.dt);
figure;plot(highpassdata')
%%
spikethresh = 4;
do_medfilt =1;
%% plots


numcol = 2;

%use first stim as template for time - assumes all same length.
    [sigon,sigoff] = GetSigTimes(filtexpt,filtexpt.stimcond,1);
    sigdur = size(filtexpt.stimcond(1).wavs,1)/44100;
 
% get cutdata with spikes removed
cutdata = [];
for istim = 1:size(filtexpt.stimcond,2)
    stimexpt = filtesweeps(filtexpt,0,'wavnames',filtexpt.stimcond(istim).wavnames);
    sigdata = stimexpt.wc.data * 1000;
    cutdata{istim} = removeSpikes(sigdata,spikethresh,stimexpt.wc.dt);
end

    
% make raster    
figure; hold on
for istim = 1:size(expt.stimcond,2)
    stimexpt = filtesweeps(filtexpt,0,'wavnames',expt.stimcond(istim).wavnames);
    sigdata = stimexpt.wc.data * 1000;
    highpassdata=HighpassGeneral(sigdata,1/expt.wc.dt);
    [spikesmat, gausstosmooth]=getspikesmat(highpassdata,spikethresh,expt.wc.dt);
    spkvec = zeros(size(spikesmat,1),size(spikesmat,2));
    spkvec(find(spikesmat)) = 1;
    %     spkrespfig = figure;

    subplot(ceil(size(expt.stimcond,2)/numcol),numcol,istim)
    hold on
    for itrial = 1:size(spkvec,1)
        thesespks = find(spkvec(itrial,:))*expt.wc.dt;
        for ispike = 1:size(thesespks,2)
            line([thesespks(ispike),thesespks(ispike)],...
                [1*itrial,(1*itrial)+1],'color','k','LineWidth',3)
        end
    end
    set(gca,'XLim',[0,(sigoff*expt.wc.dt)+2],'YLim',[0,size(spkvec,1)+1],'XTick',[],'YTick',[])
    ylims = get(gca,'YLim');
    line([sigon*expt.wc.dt,sigon*expt.wc.dt],ylims,'color','r','LineWidth',3)
    line([sigoff*expt.wc.dt,sigoff*expt.wc.dt],ylims,'color','r','LineWidth',3)
    
    if istim==1
        title(num2str(unique(filtexpt.sweeps.Vm)))

    end
end
 set(gcf,'Position',[96         577        1772         425])

if do_medfilt ==0
% plot subthresh with spikes removed
figure; hold on
for istim = 1:size(filtexpt.stimcond,2)
    thisdata = cutdata{istim};
%     thisdata = thisdata - repmat(median(thisdata(:,1:sigon),2),1,size(thisdata,2));

    subplot(ceil(size(stimexpt.stimcond,2)/numcol),numcol,istim)
    hold on
    xtime = [1:size(thisdata,2)]*stimexpt.wc.dt;
    line(xtime,thisdata','color',[0.5 0.5 0.5])
    line(xtime,mean(thisdata)','color','k','LineWidth',3)
    axis tight
    ylims = get(gca,'YLim');
    line([sigon*expt.wc.dt,sigon*expt.wc.dt],ylims,'color','r','LineWidth',3)
    line([sigoff*expt.wc.dt,sigoff*expt.wc.dt],ylims,'color','r','LineWidth',3)
    ylabel(stimexpt.stimcond(istim).wavnames)
    set(gca,'XTick',[])
    
    if istim==1
        title(num2str(unique(filtexpt.sweeps.Vm)))

    end
end
set(gcf,'Position',[96          27        1772         556])
end

if do_medfilt ==1
% plot subthresh with spikes removed and median Vm removed
figure; hold on
for istim = 1:size(filtexpt.stimcond,2)
    thisdata = cutdata{istim};
    thisdata = thisdata - repmat(median(thisdata(:,1:sigon),2),1,size(thisdata,2));


    subplot(ceil(size(stimexpt.stimcond,2)/numcol),numcol,istim)
    hold on
    xtime = [1:size(thisdata,2)]*stimexpt.wc.dt;
        if plot_stimperiod ==1
    line(xtime(sigon:sigoff),thisdata(:,sigon:sigoff)','color',[0.5 0.5 0.5])
    line(xtime(sigon:sigoff),mean(thisdata(:,sigon:sigoff))','color','k','LineWidth',3)
        end
        if plot_stimperiod ==0
    line(xtime,thisdata','color',[0.5 0.5 0.5])
    line(xtime,mean(thisdata)','color','k','LineWidth',3)
        end
    axis tight
    ylims = get(gca,'YLim');
    line([sigon*expt.wc.dt,sigon*expt.wc.dt],ylims,'color','r','LineWidth',3)
    line([sigoff*expt.wc.dt,sigoff*expt.wc.dt],ylims,'color','r','LineWidth',3)
    ylabel(stimexpt.stimcond(istim).wavnames)
    set(gca,'XTick',[])
    if istim==1
        title(num2str(unique(filtexpt.sweeps.Vm)))

    end
end
set(gcf,'Position',[96          27        1772         556])
end