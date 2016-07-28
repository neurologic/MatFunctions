cd('/Users/kperks/mnt/cube/Ice/kperks/B1087/MetaMat/')
% currentpath = pwd

% path_to_phy = '../klusta/phy041316/';

load('MVRegressPrep.mat','SignalFiltR','include_inds','isNarrow');%,'ExtracellStimNames','IntracellStimNames')

% CellID = 1; %which intracellular cell

% load('IntracellularData.mat','IntracellularData')

% load('metatoes.mat','allchanxz','metatoes')
% load('behavior.mat')
% loc = allchanxz(include_inds);


NarrowInds = isNarrow(include_inds);
Pop = SignalFiltR;

calcstart = round(size(Pop,3)/6);
fs = 1/IntracellularData{1}.dt;

Narrow = -squeeze(mean(Pop(find(NarrowInds==1),istim,:),1));
Wide = squeeze(mean(Pop(find(NarrowInds==0),istim,:),1));

figure;
hold on;

norm = max(max(Narrow));%,1,size(Narrow,2));
offset = min(min(Narrow));%,1,size(Narrow,2));
% line([1:size(Pop,3)]/fs,(Narrow./norm)-offset,'color','r')
line([1:calcstart]/fs,Narrow(:,1:calcstart),'color','r')
norm = max(max(Wide))
offset = min(min(Wide));
% line([1:size(Pop,3)]/fs,(Wide./norm)-offset,'color','b')
line([1:calcstart]/fs,Wide(:,1:calcstart),'color','b')

figure;
hold on
line([1:calcstart]/fs,mean(Narrow(:,1:calcstart)),'color','r')
line([1:calcstart]/fs,mean(Wide(:,1:calcstart)),'color','b')

figure;
hold on
line([1:size(Pop,3)]/fs,(mean(Narrow)-min(mean(Narrow)))./max(mean(Narrow)-min(mean(Narrow))),'color','r')
line([1:size(Pop,3)]/fs,(mean(Wide)-min(mean(Wide)))./max(mean(Wide)-min(mean(Wide))),'color','b')
% 
figure;
hold on
line([1:size(Pop,3)]/fs,mean(Narrow)./min(mean(Narrow)),'color','r')
line([1:size(Pop,3)]/fs,mean(Wide)./max(mean(Wide)),'color','b')

figure;
hold on
line([1:size(Pop,3)]/fs,(Narrow)-max((Narrow)),'color','r','LineWidth',2)
line([1:size(Pop,3)]/fs,(Wide)-(mean(Wide)),'color','b','LineWidth',2)
set(gca,'TickDir','out','YTick',[])
axis tight


[r,lags] = xcorr(mean(Narrow),mean(Wide));
figure;line(lags,r);
lagtime = lags(find(r == max(r)))/fs