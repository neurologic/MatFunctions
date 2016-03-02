function IntracellularData_entry = PopulateIntracellularData(filtexpt,spikethresh)

[sigon,sigoff] = GetSigTimes(filtexpt,filtexpt.stimcond,1);
sigdur = size(filtexpt.stimcond(1).wavs,1)/44100;

cutdata = [];
prestim = [];
stim = [];

for istim = 1:size(filtexpt.stimcond,2)
    stimexpt = filtesweeps(filtexpt,0,'wavnames',filtexpt.stimcond(istim).wavnames);
    sigdata = stimexpt.wc.data * 1000;
    thisdata = removeSpikes(sigdata,spikethresh,stimexpt.wc.dt);
    thisdata = thisdata - repmat(median(thisdata(:,1:sigon),2),1,size(thisdata,2));
    cutdata{istim} = thisdata;
    
    prestim{istim} = cutdata{istim}(:,1:sigon);
    stim{istim} = cutdata{istim}(:,sigon:sigoff);
    poststim{istim} = cutdata{istim}(:,sigoff:end);
end

IntracellularData_entry.exptname = filtexpt.name;    %  = expt.name
IntracellularData_entry.IC = unique(filtexpt.sweeps.Vm); %holding current
IntracellularData_entry.dt = filtexpt.wc.dt;    %  expt.wc.dt
stimcondCell = squeeze(struct2cell(filtexpt.stimcond));
IntracellularData_entry.stimnames = stimcondCell(2,:); %name of each stimulus in order of cutdata structs
IntracellularData_entry.stimtimes = [sigon,sigoff];
IntracellularData_entry.prestim_data = prestim;    %  cutdata with {istim} entries of (trial,1:sigon samples)
IntracellularData_entry.stim_data  = stim;    % cutdata with {istim} entries of (trial,sigon:sigoff samples)
IntracellularData_entry.poststim_data = poststim;

end
