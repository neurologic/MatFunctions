%%
r = rigdef('mac');
%% import sweeps
clear expt;
expt.name = 'KP_B1035_150423_p1c1';%'KP_B709_110602_p1c5' 
exptdir = 'Volumes/kperks/NetworkCorrManipExpt/EXPT_IGOR';%[r.Dir.IgorExpt]
expt = addIgorDataNew(expt,exptdir,1);
expt.wc.data=double(expt.wc.data);
%%
save(fullfile(r.Dir.Expt,[expt.name]),'expt')

% now addIgorDataNew takes care of doing stimcond and params from data
% stored during experiment in IGOR
%% condition stimulus for getting dB info and stim length as part of sweeps info
    
%%
 [expt,table,hfig]=QueryExpt(expt,r);


%%

AfterMotDur=1.5;
sweepstop = round((params.waveonset_time+AfterMotDur)/expt.wc.dt)+round(size(stimcond(1).wavs,1)/44100/expt.wc.dt);
 expt.wc.data=expt.wc.data(:,1:sweepstop);

%% 
stimcond=Acquisition_Importing_Signals(expt);
expt.stimcond=stimcond;
%% define parameters and conditions
% this could go into addIgorDataNew as well
params.steptime=round([0.05 0.3]/expt.wc.dt);
params.waveonset_time=[1.8];
params.baselinewin=round([0.35 1.8]/expt.wc.dt);
params.testAmp=-10;
params.dt=expt.wc.dt;
expt.analysis.params=params;
%% filter out any trials that were really abberant
expt=filtesweeps(expt,0,'trial',[1:32]);
      


%% save Expt
expt.wc.data=double(expt.wc.data);
save(fullfile(r.Dir.Expt,expt.name),'expt')


%% plot all trials at each holding

Vmexpt=filtesweeps(expt,0,'Vm',0);

hfig=plotallVm(Vmexpt,expt.analysis.params.steptime,expt.analysis.params.testAmp);

%% need to keep song stims with expt now even though bigger
% % need to make a stimcond for each clamp... then can use later to check if
% % same
% 
% % should really be pulling stimuli in from igor....
% % then in stimcond can check/record which holding were used for those
% % stimuli
% 
%     Vm=999;
%     sortstims.fieldnames='Vm';
%     sortstims.sortvalues=Vm;
%     folder=[r.Dir.IgorExpt expt.name(1:end-5)];
%     nrepsstim=1;
%     stimcond = filtstimcond(expt,folder,expt.analysis.params.waveonset_time,nrepsstim,sortstims);
% 
% end
% expt.stimcond=stimcond;