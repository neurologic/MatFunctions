function hfig=VisualizeForQuery(expt, saveme,r)
% r=rigdef('manu');
hfig=[];
    for ivm=1:max(size(expt.wc.allVm))
        thisdata=filtesweeps(expt,0,'Vm',expt.wc.allVm(ivm));
        
        hfig(ivm)=figure;
        hold on;
        set(hfig,'Position',[0,0,1200,650])
        %  include in figure:
        
        %     plot of data
        s1=subplot(1,2,1);
        hold on
        
        line([1:size(thisdata.wc.data,2)]*expt.wc.dt,thisdata.wc.data');
        line([1:size(thisdata.wc.data,2)]*expt.wc.dt,mean(thisdata.wc.data)','color','k','LineWidth',3);
        axis tight
        ylims=get(gca,'YLim');
        sigtime=max(size(expt.stimcond(1).wavs))/44100;
        % need a way here to account for (or at least check if) signals were all
        % the same length
        SigTimeBox(gca, expt.analysis.params.waveonset_time, ...
            (expt.analysis.params.waveonset_time+sigtime), ylims,'k')
        
        % calculate Rseries and Rinput and include on graph as well
        
        %     clamp
        title([expt.name ' clamp=' num2str(expt.wc.allVm(ivm))],'Interpreter','none');
        
        %     how many signals played in this condition?
        sigs=unique(expt.sweeps.wavnames);
        
        %     how many trials of each signal?
        trialpersig=[];
        for isig=1:max(size(sigs))
            sigexpt=filtesweeps(thisdata,0,'wavnames',sigs(isig));
            trialpersig(isig)= size(sigexpt.wc.data,1);
        end
        
        cnames={'trials'};
        rnames=sigs;
        hold on
        s2=subplot(1,2,2);
        hold on
        set(s2,'Position',[0.7, 0.1, 0.3, 0.85])
        set(s1,'Position',[0.05,0.1,0.6, 0.85])
        set(s2, 'XTick',[],'YTick',[]);
        annotation('textbox',[0.7, 0.1, 0.15, 0.85],'string',rnames);
        annotation('textbox',[0.9, 0.1, 0.1, 0.85],'string', num2str(trialpersig'));
        % t=uitable('Parent',hfig,'Data',trialpersig','ColumnName',cnames,'RowName',rnames);
        %  set(t,'Position',[800,25,300,600])
        
        %    save this figure in the folder to reference later if need to
        if saveme==1
        saveas(hfig(ivm),[r.Dir.Expt expt.name '\allTrials' num2str(expt.wc.allVm(ivm)) '.fig']);
        saveas(hfig(ivm),[r.Dir.Expt expt.name '\allTrials' num2str(expt.wc.allVm(ivm)) '.jpg']);
        end
    end