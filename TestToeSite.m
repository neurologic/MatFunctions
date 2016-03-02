function TestToeSite(toedata)
% takes toedata
% >> psth plot by depth
% >> distribution of good, mua, unsorted clusters
% >> spike statistics: shape, number, isi, violations - distribution of these measures per type of cluster

% set default plot params
set(0,'DefaultAxesFontSize',20)
set(0,'DefaultTextFontSize',20)
set(0,'DefaultTextFontName','Helvetica')
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesColorOrder',[0 0 0;1 0 1;0 1 1; 0 1 0])


%%%%%%%%%%% >> psth plot by depth
allpsth = [];
tau = 0.01;
nsamps =1000;

for isite = 1:size(metatoes,1)
    
    thissite = metatoes{isite};
    
    for istim = 1:size(metatoes{isite}.stims,1);
        thisstim = thissite.stims{istim};
        
        for itrial = 1:thisstim.ntrials
            
            thisstim_toes = thisstim.toes;
            
            thistrial_toes = thisstim_toes{itrial};
            
            filteredspikes = zeros(1,nsamps);
            
            if ~isempty(thistrial_toes)
                %                 fun = filtered_response(thistrial_toes, tau);
                fun = filtered_response_synaptic(thistrial_toes, tau);
                filteredspikes = fun(linspace(-2,8.5,nsamps));
            end
            
            allpsth(isite,istim,itrial,:) = filteredspikes;
        end
    end
end

% normalize max cluster response for each cluster to 1
for icluster = 1:size(allpsth,1)
    thiscluster = allpsth(icluster,:,:,:);
    thiscluster = thiscluster./max(thiscluster(:));
    normpsth(icluster,:,:,:) = thiscluster;
end


% trial-averaged responses by depth
maxplot_ht = 200;
scaledpsth = normpsth.*maxplot_ht;

include_inds = find(allsortclass);

% include_inds = find(allsortclass == 3);

for istim = [1,10,20,30,42:44];
    figure;
    set(gcf,'Position',[680 83 449 1015])
    hold on
    
    xtime = linspace(-2,8.5,nsamps);
    for isite = 1:size(include_inds,1)
        
        dcoff = allchanz(include_inds(isite));
        
        meanpsth = mean(squeeze(scaledpsth(include_inds(isite),istim,:,:))) + dcoff;
        
        line(xtime,meanpsth')
        
    end
    line(xtime,(50*(mean(squeeze(mean((squeeze(scaledpsth(include_inds,istim,:,:))),2)),1))+(site4_z0+200))','LineWidth',1)
    axis tight
    
    fs = 31250;
    stimfs = 40000;
    stimdur = (metatoes{isite}.stims{istim}.stim_end_times(1)/fs)-(metatoes{isite}.stims{istim}.stim_start_times(1)/fs);
    line([0,0],[site1_z0-800,site4_z0])
    line([stimdur/6,stimdur/6],[site1_z0-800,site4_z0])
    line([stimdur,stimdur],[site1_z0-800,site4_z0])
end

% for all stim, plot psth for all trials in black and trial-averaged psth in red
% nstim = size(toedata{1}.stims,1);
% ncol = 2;
% figure;
% hold on
% for istim = 1:nstim
%     % istim = 16;
%     
%     subplot(ceil(nstim/ncol), ncol, istim)
%     meanpsth_across_clusters = [];
%     for itrial = 1:20
%         %     for istim = 1:16
%         meanpsth_across_trials = squeeze(toedata_singleSite_psth(include_inds,istim,itrial,:));
%         meanpsth_across_clusters(istim,itrial,:) = mean(meanpsth_across_trials,1);
%         line(xtime,mean(meanpsth_across_trials,1),'color','k');
%         %     end
%     end
%     line(xtime,squeeze(mean(meanpsth_across_clusters,2)),'color','r')
%     axis tight
%     
%     if istim ~= nstim
%         set(gca,'XTickLabel',[])
%     end
% end

end