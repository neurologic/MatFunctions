function stimcond=Acquisition_Importing_Signals(expt);
% modified from addIgorDataNew which brings in expt struct
% this function will find the signal files from this experimentfolder make
% a signal struct instead of using filtstimcond function to make stimcond
% with wavs from expt file
r = rigdef('z');
ind = regexp(expt.name,'_');
exptfolder = expt.name(1:ind(end)-1);
exptfilename = [expt.name];

d = dir([fullfile(r.Dir.IgorExpt,exptfolder,exptfilename) '*sigs*']);
index=1;
for ifile = 1:size(d,1)
    disp(['Loading ' fullfile(r.Dir.IgorExpt,exptfolder,d(ifile).name)])
    temp = importdata(fullfile(r.Dir.IgorExpt,exptfolder,d(ind_datafiles(ifile)).name));
    nsweeps_infile = temp.data(4,end);
    firstsw = temp.data(5,end)+1;
    endsig=temp.data(6,end)*temp.data(3,end);
    for isig=1:nsweeps_infile
        stimcond(index).wavs=temp.data(1:endsig,isig); % last column is sample sweep info
        stimcond(index).wavnames=temp.textdata{1,isig}(2:end-1);
        index=index+1;
    end
    
end
