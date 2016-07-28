
cd('/Users/kperks/mnt/cube/Ice/kperks/B1083/MetaMat')
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
%% plot norm mean cluster response by location
these_classinds = find(allsortclass == 2);
these_siteinds = find(allsiteID);%find(allsiteID == 1);
include_inds = intersect(these_classinds,these_siteinds);

ntrials = size(allSpikes_SynFilter,3);
nsamps = size(xtime_spikes,2);
for istim = 1%:size(metatoes{1}.stims,1)
    hfig = figure;
    set(gcf,'Position',[680 83 449 1015],'Visible','off')
    hold on
    
%     xtime = linspace(-2,8.5,nsamps);
    for icluster = 1:size(include_inds,1)
        
        dcoff = allchanxz(include_inds(icluster));
        
        ClusterFun = squeeze(allSpikes_SynFilter(include_inds(icluster),istim,:));
        
        ClusterFiltR = zeros(ntrials,nsamps);
        for itrial = 1:ntrials
            trialfun = ClusterFun{itrial};
            if ~isempty(trialfun)
                ClusterFiltR(itrial,:) = trialfun(xtime_spikes);
            end
        end
        meanpsth = mean(ClusterFiltR,1);
        meanpsth = ((meanpsth./max(meanpsth))*50) - dcoff;
        
        line(xtime_spikes,meanpsth')
        
    end
    
    ylims = get(gca,'YLim');
    line([xtime_spikes(sigon),xtime_spikes(sigon)],ylims,'color','r')
    line([xtime_spikes(sigoff),xtime_spikes(sigoff)],ylims,'color','r')
    set(gcf,'Visible','on')
    
end

saveas(hfig,'ClusterTrialAvg_ByLocation.jpg')
saveas(hfig,'ClusterTrialAvg_ByLocation.eps')
%% plot norm mean pop response by location

these_classinds = find(allsortclass == 2);
these_siteinds = find(allsiteID == 1);
include_inds = intersect(these_classinds,these_siteinds);

ntrials = size(allSpikes_SynFilter,3);
nsamps = size(xtime_spikes,2);

locbins = unique(allchanxz(include_inds));
for istim = 2%:size(metatoes{1}.stims,1)
    figure;
    set(gcf,'Position',[680 83 449 1015],'Visible','off')
    hold on
    
    for iloc = 1:size(locbins,1)
        %     xtime = linspace(-2,8.5,nsamps);
        locpsth = [];
        loc_inds = intersect(include_inds,find(allchanxz == locbins(iloc)));
        for icluster = 1:size(loc_inds,1)
            
            dcoff = allchanxz(loc_inds(icluster));
            
            ClusterFun = squeeze(allSpikes_SynFilter(loc_inds(icluster),istim,:));
            
            ClusterFiltR = zeros(ntrials,nsamps);
            for itrial = 1:ntrials
                trialfun = ClusterFun{itrial};
                if ~isempty(trialfun)
                    ClusterFiltR(itrial,:) = trialfun(xtime_spikes);
                end
            end
            meanpsth = mean(ClusterFiltR,1);
            meanpsth = ((meanpsth./max(meanpsth))*50) + dcoff;
            
            locpsth(icluster,:) = meanpsth;
            
            
        end
        line(xtime_spikes,mean(locpsth)')
    end
    set(gcf,'Visible','on')
    
end