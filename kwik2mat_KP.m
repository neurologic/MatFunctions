
function toedata = kwik2mat_KP
%run from the site folder and it will get cluster ID and SIteID from pwd
% needs to be run from site folder to get indx_port_site.txt
% execute /Applications/MATLAB_R2015a.app/bin/matlab -nodesktop
clear toedata

% kwik to mat.  Get good units out of sorted kwik file with stimulus
% information and store in one well formatted mat file
% Brad Theilman August/September 2015
% modified KP Jan2016

% klustaID='phy012716';
% SiteID='Pen01_Lft_AP200_ML800__Site01_Z2200__B1235_cat_P01_S01_1';
currentpath = pwd;
slashinds = regexp(currentpath,'/');
klustaID = currentpath(slashinds(9)+1:slashinds(10)-1);
SiteID = currentpath(slashinds(10)+1:end);
scoreinds = regexp(SiteID,'_');
kwikname=[SiteID(scoreinds(8)+1:end) '.kwik'];
BirdID = SiteID((scoreinds(8)+1):(scoreinds(9)-1));
midnight_hit= 0;
% 


datapath = '/Users/kperks/mnt/cube/Ice/kperks/';
% datapath='/Volumes/BookDrive/Github/NetworkCorrManipExpt/Ice/rsync/';
kwikfile=[datapath BirdID '/klusta/' klustaID '/' SiteID '/' kwikname];
outfile=[datapath BirdID '/klusta/' klustaID '/' SiteID '/toedata_' SiteID '.mat'];

% h5disp(kwikfile)
% kwikinfo = h5info(kwikfile);

% Get cluster and time information for each spike
%%%%%%%%%% KNOWN ISSUE
%%%%%%%%%% that if have done any merging, there will be
%%%%%%%%%% clusters without any events assigned to it, with this format,
%%%%%%%%%% these clsuters are being added to toedeata... if go the other
%%%%%%%%%% direction of events and then find assigned cluster, would not
%%%%%%%%%% end up with empty clusters in toedata
spike_clusters = h5read(kwikfile, '/channel_groups/0/spikes/clusters/main');
spike_times = h5read(kwikfile, '/channel_groups/0/spikes/time_samples');

% extract all the unique cluster ids:
clusters = unique(spike_clusters);

% get all cluster info (1:MUA,2:good,3:unsorted,0:noise
cluster_classes = zeros(1, length(clusters));
for i = 1:length(clusters)
    
    cluster_num = clusters(i);
    cluster_attr_path = strcat('/channel_groups/0/clusters/main/', num2str(cluster_num));
    cluster_group = h5readatt(kwikfile, cluster_attr_path, 'cluster_group');
    cluster_classes(i) = cluster_group;
    
end

% Now, pull spikes from all clusters
all_spikes = cell(length(clusters), 2);
for i = 1:length(clusters)
    cluster_num = clusters(i);
    %pull spikes from teh good cluster
    this_cluster_spiketimes = spike_times(spike_clusters == cluster_num);
    all_spikes{i, 1} = cluster_num;
    all_spikes{i, 2} = this_cluster_spiketimes;
end


% Now get stimulus information
digmark_timesamples = h5read(kwikfile, '/event_types/DigMark/time_samples');
digmark_codes = cell2mat(h5read(kwikfile, '/event_types/DigMark/codes'));
stim_timesamples = h5read(kwikfile, '/event_types/Stimulus/time_samples');
stim_codes = h5read(kwikfile, '/event_types/Stimulus/codes');
stim_text = h5read(kwikfile, '/event_types/Stimulus/text');

stim_start_times = double(digmark_timesamples(digmark_codes == '<'));
stim_end_times = [double(digmark_timesamples(digmark_codes == '>'))];
intertrial_start_times = [double(digmark_timesamples(digmark_codes == '('))];
intertrial_end_times = [double(digmark_timesamples(digmark_codes == ')'))];

if size(stim_text,1)/2 ~= size(stim_start_times,1)
    midnight_hit =1;
end

if midnight_hit ==1
    cell_date_inds = regexp(stim_text,'date_');
    findinds = cellfun(@isempty,cell_date_inds,'UniformOutput',true);
    boul = ~findinds;
    date_ind = find(boul); %this is the datestamp that is messing up the stimid flow
    stim_text = [stim_text(1:date_ind-1); stim_text(date_ind+1:end)];
    stim_codes = [stim_codes(1:date_ind-1); stim_codes(date_ind+1:end)];
    stim_timesamples = [stim_timesamples(1:date_ind-1); stim_timesamples(date_ind+1:end)];
    % even = stim_text(2:2:length(stim_text));
    % stim_text = even;
end

stim_start_filename = cell2mat(stim_text);
stim_start_filename = stim_start_filename(2:2:end, :);
stim_start_filename = cellstr(stim_start_filename);

stim_files_unique = unique(stim_start_filename);
nstims = length(stim_files_unique);

%
% For each stimulus, get number of trials for that stimulus
% start building final_mat
final_data = struct;
final_data.nstims = nstims;
numtrials = zeros(1, nstims);
stim_data = cell(nstims, 1);

for i = 1:nstims
    
    stim_entry = struct;
    numtrials(i) = sum(strcmp(stim_files_unique(i), stim_start_filename));
    
    stim_start_times_this_stim = stim_start_times(strcmp(stim_files_unique(i), stim_start_filename));
    stim_end_times_this_stim = stim_end_times(strcmp(stim_files_unique(i), stim_start_filename));
    stim_entry.name = stim_files_unique(i);
    stim_entry.start_times = stim_start_times_this_stim;
    stim_entry.end_times = stim_end_times_this_stim;
    stim_entry.ntrials = sum(strcmp(stim_files_unique(i), stim_start_filename));
    stim_data{i, 1} = stim_entry;
    
end

% Go through each cell and build final data matrix

% fs = 31250.0; 
% get fs from data file itself 
a = h5info(kwikfile,'/application_data/spikedetekt/');
a = struct2cell(a.Attributes);
attribute_names = a(1,:);
ind = find(strcmp(attribute_names,'sample_rate'));
fs = a{4,ind};

pre_stim_duration = 2; %in seconds
post_stim_duration = 2; %seconds
pre_stim_duration_samps = pre_stim_duration*fs;
post_stim_duration_samps = post_stim_duration*fs;

n_units = length(clusters);
toedata = cell(n_units, 1);



% load in cluster probe location information
clusta_info_file = [datapath BirdID '/klusta/' klustaID '/' SiteID '/cluster_info.mat'];
clusta_info_mat = load(clusta_info_file);


for unit_num = 1:length(clusters)
    
    unit_entry = struct;
    unit_entry.id = clusters(unit_num);
    unit_entry.sort_class = cluster_classes(unit_num);
    
    primary_location_prb_microns = clusta_info_mat.prb_location_microns(find(clusta_info_mat.cluster_id == clusters(unit_num)),:);
    
    unit_entry.primary_location_prb_microns = primary_location_prb_microns;
    unit_entry.fs = fs;
    stims = cell(nstims, 1);
    
    spiketimes_thisunit = double(all_spikes{unit_num, 2});
    % Loop through each stimulus
    for stim_num = 1:nstims
        
        stim_unit_entry = struct;
        
        this_stim_data = stim_data{stim_num, 1};
        stim_unit_entry.name = this_stim_data.name;
        stim_unit_entry.ntrials = this_stim_data.ntrials;
        
        %create trial times
        trial_start_times = this_stim_data.start_times - pre_stim_duration_samps;
        trial_end_times = this_stim_data.end_times + post_stim_duration_samps;
        stim_unit_entry.stim_start_times = this_stim_data.start_times;
        stim_unit_entry.stim_end_times = this_stim_data.end_times;
        stim_unit_entry.trial_start_times = trial_start_times;
        stim_unit_entry.trial_end_times = trial_end_times;
        
        %Go through each trial to divy up spike times relative to stim
        %start
        toes = cell(this_stim_data.ntrials, 1);
        for trialnum = 1:this_stim_data.ntrials
            trial_start = trial_start_times(trialnum);
            trial_end = trial_end_times(trialnum);
            
            
            
            spiketimes_samps_thistrial = spiketimes_thisunit(spiketimes_thisunit >= trial_start & spiketimes_thisunit <= trial_end);
            spiketimes_samps_thistrial_relstimonset = spiketimes_samps_thistrial - this_stim_data.start_times(trialnum);
            spiketimes_secs_thistrial_relstimonset = double(spiketimes_samps_thistrial_relstimonset) / fs;
            
            toes{trialnum, 1} = spiketimes_secs_thistrial_relstimonset;
        end
    
        stim_unit_entry.toes = toes;
        stims{stim_num, 1} = stim_unit_entry;     
    end
    
    unit_entry.stims = stims;
    unit_entry.all_spikes = double(all_spikes{unit_num, 2});
    toedata{unit_num, 1} = unit_entry;
    
   
end


% Format output file name

% data_to_save = struct();
% data_to_save.birdID = '';
% data_to_save.penetrationID = 1;
% data_to_save.siteID = 3;
% data_to_save.XYpos = [0, 0];
% data_to_save.Zpos = 0;
% data_to_save.target_structure = 'NCM';
% 
% data_to_save.fs = fs;
% data_to_save.toedata = toedata;
% 
% sav_date = datestr(now, 30);

save(outfile, 'toedata');

% end

