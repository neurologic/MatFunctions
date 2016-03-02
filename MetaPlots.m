
%% set default plot params
set(0,'DefaultAxesFontSize',20)
set(0,'DefaultTextFontSize',20)
set(0,'DefaultTextFontName','Helvetica')
set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesColorOrder',[0 0 0;1 0 1;0 1 1; 0 1 0])
%%

sortID = 2;

nstims = size(metatoes{1}.stims,1);

hfig = figure('Visible','off');
% figname = ['Stim_' toedata{1}.stims{istim}.name{:}];
set(hfig,'Position',[ 176          64        1611        1028])
ncol = 4;


for istim = 1:nstims
    % Get the data for the chose cell
    stimnum = istim;
    subplot(ceil(nstims/ncol), ncol, istim)
    line_ind = 0;
    for isite = 1:size(metatoes,1)
        if metatoes{isite}.sort_class == 1 || metatoes{isite}.sort_class == 2
            ys = [0 + line_ind, 1+line_ind];
            unit_data = metatoes{isite, 1};
            stim_data = unit_data.stims{stimnum, 1};
            stim_end_secs = double(stim_data.stim_end_times - stim_data.stim_start_times)/fs;
            ntrials = stim_data.ntrials;
            
            trialnum = 1;%:ntrials
            
            if ~isempty(stim_data.toes{trialnum, 1})
                for spikenum = 1:length( stim_data.toes{trialnum, 1})
                    line([stim_data.toes{trialnum, 1}(spikenum), stim_data.toes{trialnum, 1}(spikenum)], ys);
                end
            end
            
            line_ind = line_ind+1;
        end
        
    end
    xlim([-2, stim_end_secs(trialnum)+2])
    ylim([0, line_ind+1]);
    line([0, 0], [0, line_ind+1], 'Color', 'red');
    line([stim_end_secs(trialnum), stim_end_secs(trialnum)], [0, line_ind+1], 'Color', 'red');
    set(gca,'TickDir','out','YTick',[1,line_ind-1],'YTickLabel',[1,line_ind-1],'XTick',[-2:2:8],'XTickLabel',[-2:2:8])
    ylabel('UnitID')
    xlabel('seconds')
    %     saveas(hfig,[outfile(1:end-11) figname],'png')
    
    %     close(hfig)
    
end
set(hfig,'Visible','on')
