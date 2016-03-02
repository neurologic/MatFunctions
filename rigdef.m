function r = rigdef(location)
if regexp(location,'usbhome')
    r.Dir.Expt = 'E:\data\experiments\';
    r.Dir.Wavs = 'C:\krista\KPstims\singlemotiflibrary_44100Hz\';
end
if regexp(location,'manu')
    r.Dir.IgorExpt = 'C:\KPexperiments\';
    r.Dir.Expt = 'C:\krista\data\experiments\';
    r.Dir.Wavs = 'C:\krista\KPstims\AllStimToPlot\';
    r.Dir.Behavior = 'C:\krista\data\behavior\';
    r.Dir.Calibrate = 'C:\krista\KPstims\AllStimToPlot\Calibrate\';
end
if regexp(location,'z')
    %this is when lintu is a mapped network drive and i am working under my
    %('krista') AD account
    r.Dir.IgorExpt = 'Z:\KPexperiments\';
    r.Dir.Expt = 'Z:\krista\data\experiments\';
    r.Dir.Wavs = 'Z:\krista\KPstims\AllStimToPlot\';
    r.Dir.Behavior = 'Z:\krista\data\behavior\';
    r.Dir.Calibrate = 'Z:\krista\KPstims\Calibrated\';
end
if regexp(location,'mac')
    r.Dir.Expt = '/Users/kperks/GitHub/Data_Mat/';
%     r.Dir.IgorExpt = 'smb://lintu.ucsd.edu/kperks/KPexperiments/';
    r.Dir.EJNmsRepo_expt = '/Users/kperks/GitHub/PerksGentner2014/WholeCellData_repexpts/';
    r.Dir.IgorExpt = '/Volumes/kperks/KPexperiments/';
    r.Dir.Stims = '/Users/kperks/GitHub/StimsAllWavs/';
end

end