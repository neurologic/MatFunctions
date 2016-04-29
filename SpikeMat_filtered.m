function allSpikes = SpikeMat_filtered(toedata,type)
% type is 'gaussian' or 'synaptic'
% edges = [trialstart_time,trialend_time] 
% trialstart_time = -2  because spikes from toedata are relative time to stimulus onset
% trialend_time = 8.5 because stimuli are 6s long right now for this expt

%build psth matrix for a given stimulus on a given trial
%psth are filtered by gaussian
allSpikes = cell(size(toedata,1),size(toedata{1}.stims,1),size(toedata{1}.stims{1}.toes,1));
tau = 0.01;


for isite = 1:size(toedata,1)
    
    thissite = toedata{isite};
    
    for istim = 1:max(size(toedata{isite}.stims));
        thisstim = thissite.stims{istim};
        
        for itrial = 1:thisstim.ntrials
            
            thisstim_toes = thisstim.toes;
            
            thistrial_toes = thisstim_toes{itrial};
            
%             filteredspikes = zeros(1,nsamps);
            
            if ~isempty(thistrial_toes)
                if strcmp(type,'gaussian')
                    fun = filtered_response(thistrial_toes, tau);
                    allSpikes{isite,istim,itrial,:} = fun;
                end
                if strcmp(type,'synaptic')
                    fun = filtered_response_synaptic(thistrial_toes, tau);
                    allSpikes{isite,istim,itrial,:} = fun;
                end
%                 filteredspikes = fun(linspace(trial_edges(1),trial_edges(2),nsamps));
            end
            
%             allSpikes(isite,istim,itrial,:) = filteredspikes;
        end
    end
end