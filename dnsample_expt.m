function expt=dnsample_expt(expt,newSampleRate)
olddt=expt.wc.dt;
bin=round((1/olddt)/newSampleRate);
expt.wc.dt=olddt*bin;
expt.wc.data=expt.wc.data(:,1:bin:end);
end