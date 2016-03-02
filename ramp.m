function rampedY=ramp(Y,rampms)
numrampsamp = floor(rampms/1000*44100);
ramp = ([0:numrampsamp-1]/numrampsamp)';
rampedY = [Y(1:numrampsamp).*ramp ; Y(numrampsamp+1:end-numrampsamp); Y(end-(numrampsamp-1):end).*flipud(ramp)];

