function [expt,table,hfig]=QueryExpt(expt,r)

% Vm
% which clamps do you care about for this cell
% if there are any you dont care about, then filter them out and chamge
% expt.wc.allVm to reflect that change
hfig=VisualizeForQuery(expt,0,r);

% excluded Trials
% exclude those trials
trialsexclude=input('are there any trials you want to exclude from analysis?');
expt.trialsexcluded=trialsexclude;
if ~isempty(trialsexclude)
    alltrials=expt.sweeps.trial;
    alltrials(trialsexclude)=0;
    keepinds=find(alltrials);
    keeptrials=alltrials(keepinds);
    expt=filtesweeps(expt,0,'trial',keeptrials);
end

clamps=input('what clamp values do you want to use?');
excludeclamp=[];
for iclamp=1:max(size(expt.wc.allVm))
    index=1;
    if isempty(find(clamps==expt.wc.allVm(iclamp)))
        excludeclamp(index)=expt.wc.allVm(iclamp);
        index=index+1;
    end
end
expt.excludedclamps=excludeclamp;
if ~isempty(excludeclamp)
    expt=filtesweeps(expt,0,'Vm',clamps);
end
expt.wc.allVm=unique(expt.sweeps.Vm);
for ifig=1:max(size(hfig))
    close(hfig(ifig));
end

%  re-save experiment
save(fullfile(r.Dir.Expt,expt.name),'expt')

% now for each clamp potential, replot, save, and go through and document info:
hfig=VisualizeForQuery(expt,1,r);

%  from mean response at each clamp, calculate Rin and Rs and provide a
%  plot of the step to save
for iclamp=1:max(size(clamps))
    table(iclamp).clamp=clamps(iclamp);
    vmexpt=filtesweeps(expt,0,'Vm',clamps(iclamp));
    whatclamp=unique(vmexpt.sweeps.clamp);
    assert((max(size(whatclamp))==1),'different clamp types reported for this one holding');
    isIC=whatclamp;
    [rin,rs]=WholeCellParams(vmexpt,isIC); %in Igor, when current clamp box checked, a "1" is stored as the clamp for that trial and 0 is for VC
    table(iclamp).Rin=rin;
    table(iclamp).Rs=rs;
    
    %holding current or resting potential:
    baselinedata=mean(vmexpt.wc.data(:,vmexpt.analysis.params.baselinewin(1):vmexpt.analysis.params.baselinewin(2)),1);
    holdingvalue=median(baselinedata);
    table(iclamp).baselinedata=baselinedata;
    table(iclamp).holdingvalue=holdingvalue;
    
    % type
    %  1:vme 2:vmi 3:IC 4:spikes
    typestr={'vme','vmi','ic','spikes'};
    printstr=['for clamp=' num2str(clamps(iclamp)) ' __  1:vme 2:vmi 3:IC 4:spikes   '];
    type=input(printstr);
    table(iclamp).clamptype=typestr{type};
    
    % filter
    % if type=4 filter is high, else filter is low
    if type==4
        table(iclamp).filter='h';
    else  table(iclamp).filter='l';
    end
    
    % signals played at this clamp
    table(iclamp).sigsplayed=unique(vmexpt.sweeps.wavnames);
    
    % trials per signal recorded
    numtrials=[];
    for isig=1:length(table(iclamp).sigsplayed)
        sigexpt=filtesweeps(vmexpt,0,'wavnames', table(iclamp).sigsplayed{isig});
        numtrials(isig)=size(sigexpt.wc.data,1);
    end
    table(iclamp).numtrials=numtrials;
    
end
close(hfig)
expt.table=table;
save(fullfile(r.Dir.Expt,expt.name),'expt')
