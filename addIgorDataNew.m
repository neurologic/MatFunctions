function expt = addIgorDataNew(expt,exptdir,downsample)
 % this function will find all the files from this exptmentfolder and files
 % and concatenate them
 
%  r = rigdef('mac');
 % get exptfolder
 ind = regexp(expt.name,'_');
    exptfolder = expt.name(1:ind(3)-1);
%    exptfolder = expt.name(1:ind(end-1)-1);
%    exptfolder = expt.name(1:ind(end)-1);
 exptfilename = [expt.name];
 
 % file files in experiment
 d = dir([fullfile(exptdir,exptfolder,exptfilename) '*.igor2matlab']);

 % get nsamples and sweeps etc
 for ifile = 1:length(d)
%      if ~isempty(regexp(d(ifile).name,'epochinfo'))
          if ~isempty(regexp(d(ifile).name,'epochinfo'))
         temp = importdata(fullfile(exptdir,exptfolder,d(ifile).name));
         break
     end
 end
 expt.daqinfo.nchn = temp.data(1);
 expt.daqinfo.nsamples = temp.data(2);
 expt.daqinfo.samplerate = temp.data(3);
 expt.daqinfo.nsweeps = temp.data(4);
 
notdata = ifile;

% load info about V-Clamp vs I-clamp per trial
for ifile = 1:length(d)
     if ~isempty(regexp(d(ifile).name,'clampinfo'))
         temp = importdata(fullfile(exptdir,exptfolder,d(ifile).name));
         break
     end
end
expt.sweeps.clamp = nan(expt.daqinfo.nsweeps,1);
expt.sweeps.clamp = temp.data(2:expt.daqinfo.nsweeps+1,:);
notdata = [notdata ifile];

% load info about timeStamp per trial
for ifile = 1:length(d)
     if ~isempty(regexp(d(ifile).name,'timestampinfo'))
         temp = importdata(fullfile(exptdir,exptfolder,d(ifile).name));
         break
     end
end

expt.sweeps.time = temp(3:end);
notdata = [notdata ifile];

% load info about holding pot/cur per trial
for ifile = 1:length(d)
     if ~isempty(regexp(d(ifile).name,'holdinginfo'))
         temp = importdata(fullfile(exptdir,exptfolder,d(ifile).name));
         break
     end
end
expt.sweeps.holding = nan(expt.daqinfo.nsweeps,1);
if expt.daqinfo.nsweeps > size(temp.data,1)
    expt.sweeps.holding(1:size(temp.data,1),:) = temp.data;
else
    expt.sweeps.holding = temp.data(1:expt.daqinfo.nsweeps,:); %this line was the original without the if loop
end
notdata = [notdata ifile];

% load info about Rin per trial
for ifile = 1:length(d)
     if ~isempty(regexp(d(ifile).name,'Rininfo'))
         temp = importdata(fullfile(exptdir,exptfolder,d(ifile).name));
         break
     end
end
% expt.sweeps.rin = nan(expt.daqinfo.nsweeps,1);
% expt.sweeps.rin = temp.data(1:expt.daqinfo.nsweeps,:);
notdata = [notdata ifile];

% load info about Rs per trial
for ifile = 1:length(d)
     if ~isempty(regexp(d(ifile).name,'Rsinfo'))
         temp = importdata(fullfile(exptdir,exptfolder,d(ifile).name));
         break
     end
end
% expt.sweeps.rs = nan(expt.daqinfo.nsweeps,1);
% expt.sweeps.rs = temp.data(1:expt.daqinfo.nsweeps,:);
notdata = [notdata ifile];

% load sweep Vm info
 for ifile = 1:length(d)
     if ~isempty(regexp(d(ifile).name,'vminfo'))
         temp = importdata(fullfile(exptdir,exptfolder,d(ifile).name));
         break
     end
 end
%  
% temp=111;
% expt.sweeps.Vm=repmat(temp,480,1);
 expt.sweeps.Vm = '';
 expt.sweeps.Vm= temp.data;
 %expt.sweeps.Vm= temp.data(1:expt.daqinfo.nsweeps,:);
%  a=111;
%  aa=repmat(a,10,1);
%  tempdata=[aa' temp.data']';
%  expt.sweeps.Vm= tempdata(1:expt.daqinfo.nsweeps,:);
 
 notdata = [notdata ifile];
 
% load sweep info
 for ifile = 1:length(d)
     if ~isempty(regexp(d(ifile).name,'sweepinfo'))
         temp = importdata(fullfile(exptdir,exptfolder,d(ifile).name));
         break
     end
 end
 
 %take this out!!!!!!!
%  for i=130:176
%     temp{i}='xxxxxxxxxxxxxxxxxxxxxxxxx0_00' 
%     temp{i}='009nsbEEEs001c084p001.wav0_80'
%  end
 %temp{181}=179pnsb105s010c00Xp036.wav4_00
 % add sweeps field
 %a=temp{11};
 %temp{11}=['00' a];
 
%  expt.sweeps.motifID ='';
%  expt.sweeps.wavnames ='';
 expt.sweeps.motifID =[];
 expt.sweeps.wavnames =[];
 for isweep = 3:length(temp)
     %      ind=regexpi(temp{isweep},'wav');
     %      expt.sweeps.motifID{isweep-2,1} = temp{isweep}(4:(ind+2));
     if ~isempty(regexp(temp{isweep},'ICstim'))
%          temp{isweep}=[temp{isweep} 'xxxxxxxxxx'];
            ind=regexpi(temp{isweep},'txt');
          expt.sweeps.motifID{isweep-2,1} = temp{isweep}(4:(ind+2));
     end
     if ~isempty(regexp(temp{isweep},'pnsb'))
%          temp{isweep}=[temp{isweep} 'xxxxxxxxxx'];
            ind=regexpi(temp{isweep},'wav');
          expt.sweeps.motifID{isweep-2,1} = temp{isweep}(7:(ind+2));
     end
          if ~isempty(regexp(temp{isweep},'pns'))
%          temp{isweep}=[temp{isweep} 'xxxxxxxxxx'];
            ind=regexpi(temp{isweep},'wav');
          expt.sweeps.motifID{isweep-2,1} = temp{isweep}(7:(ind+2));
     end
%      if ~isempty(regexp(temp{isweep},'testTORC'))
%          expt.sweeps.motifID(isweep-2,:) = [temp{isweep}(4:18) 'xxxxxxxx'];
%      else
%          ind=regexpi(temp{isweep},'wav');
%          expt.sweeps.motifID{isweep-2,1} = temp{isweep}(4:(ind+2));
%                   expt.sweeps.motifID{isweep-2,1} = temp{isweep}(4:26);
%      end
 end
 
 for i=1:size(expt.sweeps.motifID,1)
     ind=regexpi(temp{isweep},'wav');
     expt.sweeps.wavnames{i,1}= expt.sweeps.motifID{i}(1:end-4);
%      expt.sweeps.wavnames{i,1}= expt.sweeps.motifID(i,4:end);
 end 
 
 
 notdata = [notdata ifile];

 expt.wc.allVm=unique(expt.sweeps.Vm);
 
%add sweeps.trial numbers
expt.sweeps.trial=[1:expt.daqinfo.nsweeps]';
expt.sweeps.trialindices=[1:size(expt.sweeps.trial,1)]';

% load stim wav info
 for ifile = 1:length(d)
     if ~isempty(regexp(d(ifile).name,'wavlistinfo'))
         temp = importdata(fullfile(exptdir,exptfolder,d(ifile).name));
         break
     end
 end
 expt.song.wavlist='';
 expt.song.wavlist=temp(2:end);
 numstims=size(expt.song.wavlist,1);
notdata = [notdata ifile];


% d = dir([fullfile(exptdir,exptfolder,exptfilename) '*sigs*']);
% index=1;
for ifile = 1:size(d,1)
    if ~isempty(regexp(d(ifile).name,'sigs'))
        disp(['Loading ' fullfile(exptdir,exptfolder,d(ifile).name)])
        temp = importdata(fullfile(exptdir,exptfolder,d(ifile).name));
        nsweeps_infile = temp.data(4,end);
        firstsw = temp.data(5,end)+1;
        endsig=round(temp.data(6,end)*temp.data(3,end));
        if isnan(endsig)
            endsig=round(1.5*temp.data(3,end));
            errorstr='after motif duration set to 1.5 in Matlab';
        end
       
        for isig=1:nsweeps_infile
            stimcond(isig).wavs=temp.data(1:(end-endsig),isig); % last column is sample sweep info
            wavlength(isig)=max(size(temp.data(1:(end-endsig),isig)));
            stimcond(isig).wavnames=temp.textdata{1,isig}(5:end-5);
        end
        expt.stimcond=stimcond;
        break
    end
end
notdata = [notdata ifile];
% % load stim length info
%  for ifile = 1:length(d)
%      if ~isempty(regexp(d(ifile).name,'motiflengthinfo'))
%          temp = importdata(fullfile(exptdir,exptfolder,d(ifile).name));
%          break
%      end
%  end
%   expt.song.wavlength= temp.data(1:numstims);
% notdata = [notdata ifile];

 % load data files
 ind_datafiles = ones(size(d,1),1); %expt.daqinfo.nsweeps
 ind_datafiles(notdata) = 0;
 ind_datafiles = find(ind_datafiles);
 data = nan(expt.daqinfo.nsamples,expt.daqinfo.nsweeps,'single');
 %data = nan(132346,expt.daqinfo.nsweeps,'single');
 
 for ifile = 1:length(ind_datafiles)
     disp(['Loading ' fullfile(exptdir,exptfolder,d(ind_datafiles(ifile)).name)])
     temp = importdata(fullfile(exptdir,exptfolder,d(ind_datafiles(ifile)).name));
     nsweeps_infile = temp.data(4,end);
     firstsw = temp.data(5,end);
     
     data(:,firstsw:firstsw+nsweeps_infile-1) = temp.data(:,1:end-1); % last column is sample sweep info
 end
 expt.wc.data = data';
 expt.wc.dt=1/expt.daqinfo.samplerate;
 if downsample == 1
 expt=dnsample_expt(expt,10000);
 end
 
 
 for ifile = 1:length(d)
     if ~isempty(regexp(d(ifile).name,'epochinfo'))
         temp = importdata(fullfile(exptdir,exptfolder,d(ifile).name));
         break
     end
 end
%  EpochInfo = 
% {1-num_chans, 2-no_samples, 3-kHz*1000, 
%     4-num_sweeps,5-(bbeginstep),6-(sstepdur),
%     7-(aafterstep),8-(aaftermotif)};
 expt.analysis.params.steptime=round([temp.data(5) (temp.data(5)+temp.data(6))]/expt.wc.dt);
 expt.analysis.params.waveonset_time=temp.data(5)+temp.data(6)+temp.data(7);
 expt.analysis.params.sweepoffset_time = temp.data(8);
 %  for baselinewin... allow 20msec for RC of step to return to baseline
 expt.analysis.params.baselinewin=round([(temp.data(5)+temp.data(6)+0.02) expt.analysis.params.waveonset_time]/expt.wc.dt);
 
 %%%%%% don't know when or how this was useful, but it started to screw
 %%%%%% things up!
%  assert(max(size(unique(wavlength)))==1, 'signals of different durations');
%  AfterMotDur = expt.analysis.params.sweepoffset_time;
%  sweepstop = round((expt.analysis.params.waveonset_time+AfterMotDur)/expt.wc.dt)+round(size(expt.stimcond(1).wavs,1)/44100/expt.wc.dt);
%  expt.wc.data=expt.wc.data(:,1:sweepstop);
