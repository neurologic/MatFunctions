function [newY,meandbNew,peakdbNew,runningdb]=scalewav2(y,newdb,ismean);

NBITS=16;
bw=6.0206*NBITS;
%strip the DC offset
dcoff = (mean(y));
nodc = y-dcoff;
maxold = max(nodc);
meanrms = sqrt(mean(nodc.^2));
meandb = bw + (20*log10(meanrms));
peakrms = sqrt(max(nodc.^2));
minrms = sqrt(min(nodc.^2));
peakdb = bw + (20*log10(peakrms));
newrms = 10^((newdb-bw)/20);
mindb = bw + (20*log10(minrms));

if ismean==1
xform = 'mean';
scale = newrms/meanrms;
newY=scale*nodc;
maxnew = max(newY);

     meanrmsNew = sqrt(mean(newY.^2));
meandbNew= bw + (20*log10(meanrmsNew));
peakrmsNew = sqrt(max(newY.^2));
peakdbNew = bw + (20*log10(peakrmsNew));

runningenvelope=abs(hilbert(newY));
runningrms=sqrt(runningenvelope.^2);
runningdb=bw + (20.*log10(runningrms));

end
if ismean==0;
    xform = 'peak';
     olddb = bw + (20*log10(peakrms));
     scale = newrms/peakrms;
     newY = scale*nodc;
     maxnew = max(newY);
     
     meanrmsNew = sqrt(mean(newY.^2));
meandbNew= bw + (20*log10(meanrmsNew))
peakrmsNew = sqrt(max(newY.^2));
peakdbNew = bw + (20*log10(peakrmsNew))

runningenvelope=abs(hilbert(newY));
runningrms=sqrt(runningenvelope.^2);
runningdb=bw + (20.*log10(runningrms));
end