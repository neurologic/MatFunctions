% Krista Perks
% modified from http://klab.smpp.northwestern.edu/wiki/images/6/6b/Deneve_spike_coding_CoSMo2015.m
% http://klab.smpp.northwestern.edu/wiki/images/b/b7/DeneveCoSMo2015.pdf
% for Boerlin M. et al 2013 PLOS Computational Biology



%%%%%%%%%%%%
%ORIGINAL
%%%%%%%%%%%%
Gamma=[1,1,1,1,1,-1,-1,-1,-1,-1];
N=length(Gamma);
T=10000; Tpstart=2000; Tpend=4000; Tnstart=6000; Tnend=8000; V=zeros(N,T); O=zeros(N,T);
x=zeros(1,T); xhat=zeros(1,T); lambda=10;
s=zeros(1,T); s(Tpstart:Tpend)=200; r=zeros(N,T);
s(Tnstart:Tnend)=-200; 
epsilon=0.001;
nu=0.1;

W=Gamma'*Gamma;
NormG=diag(W)/2+nu/2;
dt=0.0001;

for t=1:T
    V(:,t+1)=(1-lambda*dt)*V(:,t)+Gamma'*s(t)*dt - W*O(:,t) + lambda*dt*W*r(:,t); 
    r(:,t+1)=r(:,t)-lambda*dt*r(:,t)+O(:,t);
    
    crit=V(:,t+1)-(NormG+epsilon*randn(N,1));
    
    O(:,t+1)=(crit>0);
    
    if (sum(O(:,t+1))>1)
        [u,v]=max(crit);
        
        O(:,t+1)=0;
        O(v,t+1)=1;
    end
    x(t+1)=x(t)+s(t)*dt;
    xhat(t+1)=(1-lambda*dt)*xhat(t)+Gamma*O(:,t);
end

plot(x','b'); hold on; plot(xhat','r'); plot(-O(2,:)+4,'r'); plot(O(1,:)+4,'g'); plot(O(3,:)+6,'m');plot(-O(4,:)+6,'k');
plot(O(1,:)+4,'g');

plot(s/5,'g')
hold off

%%
% Gamma=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1];
N=10; %number of neurons
g = [0.5:0.01:1]';
Gamma = [datasample(g,N/2);datasample(g,N/2)]'; % 1 by N vector of projection weights

[motif_data,motif_fs] = audioread('Z:\Ice\1sMotifs\Intro\01.wav');
old_t = size(motif_data,1)/motif_fs;
new_fs = 10000;
s = motif_data';
% s = hilbert(s);
% s = abs(s);
% s = medfilt1(s,200,[],2);
% s = dnsample_data(s,motif_fs,new_fs);
s = (s/max(s))*200;

T=size(s,2);  
fs_resamp = round(T/old_t);
dt = 1/fs_resamp;
% dt=0.0001;

% T = 10000;
% Tpstart=2000; Tpend=4000; Tnstart=6000; Tnend=8000;
% s=zeros(1,T); 
% s(Tpstart:Tpend)=200; 
% s(Tnstart:Tnend)=-200; 

lambdaD=10;
epsilon=0.1;
nu=0.1; 

V=zeros(N,T); 
O=zeros(N,T);
x=zeros(1,T); 
xhat=zeros(1,T); 
r=zeros(N,T);

W=Gamma'*Gamma;
NormG=diag(W)/2+nu/2;


for t=1:T
    V(:,t+1)=(1-lambdaD*dt)*V(:,t)+Gamma'*s(t)*dt - W*O(:,t) + lambdaD*dt*W*r(:,t); 
    r(:,t+1)=r(:,t)-lambdaD*dt*r(:,t)+O(:,t);
    
    crit=V(:,t+1)-(NormG+epsilon*randn(N,1));
    
    O(:,t+1)=(crit>0);  
    
    if (sum(O(:,t+1))>1)
        [u,v]=max(crit);
        
        O(:,t+1)=0;
        O(v,t+1)=1;
    end
    x(t+1)=x(t)+s(t)*dt;
    xhat(t+1)=(1-lambdaD*dt)*xhat(t)+Gamma*O(:,t);
end

figure
plot(x','b'); hold on; plot(xhat','r'); 
% plot(s/50,'g')


figure; hold on
for i = 1:N
plot(-O(i,:)+i,'k'); 
end

hfigr = figure;
hold on
hfigv = figure;
plotinds = [1,N/2,N/2+1,N];
cols = colormap(jet(size(plotinds,2)))
for i = 1:4
    figure(hfigr)
line([1:T+1]*dt, r(i,:)','color',cols(i,:))
figure(hfigv)
line([1:T+1]*dt, V(i,:)','color',cols(i,:))
end

xhat_tmp = -Gamma * r + Gamma * O;

%% 
% Gamma=[1,1,-1,-1];
% N=length(Gamma);
N=10; %number of neurons
g = 0.1; %[0.1:0.01:1]';
Gamma = [datasample(g,N/2);-datasample(g,N/2)]'; % 1 by N vector of projection weights


T=10000; 
Tpstart=2000; Tpend=4000; Tnstart=6000; Tnend=8000;
s=zeros(1,T);
s(Tpstart:Tpend)=200; 
s(Tnstart:Tnend)=-200; 

V=zeros(N,T);  
O=zeros(N,T);
r=zeros(N,T);
x=zeros(1,T);
xhat=zeros(1,T); 
lambdaD=10;
lambdaS = 100;

epsilon=0.1;
nu=0.1;

W=Gamma'*Gamma;
NormG=diag(W)/2+nu/2;
dt=0.0001;

for t=1:T
    V(:,t+1)=(1-lambdaD*dt)*V(:,t)+Gamma'*s(t)*dt - W*O(:,t) + lambdaD*dt*W*r(:,t); 
    r(:,t+1)=r(:,t)-lambdaD*dt*r(:,t)+O(:,t);
    
    crit=V(:,t+1)-(NormG+epsilon*randn(N,1));
    
    O(:,t+1)=(crit>0);
    
    if (sum(O(:,t+1))>1)
        [u,v]=max(crit);
        
        O(:,t+1)=0;
        O(v,t+1)=1;
    end
    x(t+1)=x(t)+s(t)*dt;
%     x(t+1) = (-lambdaS * x(t) * dt) + s(t)*dt + x(t);
    xhat(t+1)=(1-lambdaD*dt)*xhat(t)+Gamma*O(:,t);
end

% figure;
% hold on
% plot(x','b'); hold on; plot(xhat','r'); plot(-O(2,:)+4,'r'); plot(O(1,:)+4,'g'); plot(O(3,:)+6,'m');plot(-O(4,:)+6,'k');
% plot(O(1,:)+4,'g');
% 
% plot(s/5,'g')
% hold off
figure
plot(x','b'); hold on; plot(xhat','r'); 
plot(s,'g')


figure; hold on
for i = 1:N
plot(-O(i,:)+i,'k'); 
end

figure;
hold on
plot(r')

figure;
hold on
plot(V')


%%
load('metatoes.mat','metatoes');
load('metatoes.mat','allsortclass');
%%
goodtoes = metatoes(find(allsortclass==2));

thesetoes = goodtoes(1:20);
N=size(thesetoes,1); %number of neurons

stimstart = 0;
stimstop = (thesetoes{1}.stims{1}.stim_end_times(1) - thesetoes{1}.stims{1}.stim_start_times(1) ) / thesetoes{1}.fs;
fs = metatoes{1}.fs;
dt = 1/fs;
T = (stimstop-stimstart)*fs;

istim = 1;
xtime = round(linspace(stimstart,stimstop,(fs*(stimstop-stimstart))),4);
itrial = 1;

O = zeros(N,fs*stimstop);
for icluster = 1:N
thisstim = thesetoes{icluster}.stims{istim};
trialspikes = thisstim.toes{itrial};
stimspikes = round(trialspikes(intersect(find(trialspikes>stimstart), find(trialspikes<stimstop))),4);

O(icluster,:) = ImpulseFunction(xtime,stimspikes');
end

r=zeros(N,T);
for t=1:T-1
    r(:,t+1)=r(:,t)-lambdaD*dt*r(:,t)+O(:,t);
end
[sigma,shrinkage]=covCor(r');

sigmaR1 = sigma(1,:);
ispositive = find(sigmaR1>0);
isnegative = find(sigmaR1<0);
Gamma = diag(sigma)';
Gamma(isnegative) = -Gamma(isnegative);
% Gamma(ispositive) = -Gamma(ispositive);
% W = sigma;
W = Gamma'*Gamma;
% g = [0.9:0.01:1]';
% Gamma = [datasample(g,ceil(N/2));datasample(g,floor(N/2))]'; % 1 by N vector of projection weights

V=zeros(N,T);  
Ohat=zeros(N,T);
r=zeros(N,T);
x=zeros(1,T);
xhat=zeros(1,T); 
xhatOhat=zeros(1,T); 
shat = zeros(1,T);
lambdaD=10;

epsilon=0.01;
nu=1;

% W=Gamma'*Gamma;
NormG=diag(W)/2+nu/2;
dt=1/fs;
% dxdt = -lambda * x(t) + Gamma * o(t)
for t=1:T-1
    V(:,t+1)=(1-lambdaD*dt)*V(:,t) - W*O(:,t) + lambdaD*dt*W*r(:,t); %+Gamma'*s(t)*dt
    r(:,t+1)=r(:,t)-lambdaD*dt*r(:,t)+O(:,t);
    
    crit=V(:,t+1)-(NormG+epsilon*randn(N,1));
    
    Ohat(:,t+1)=(crit>0);
    
    if (sum(Ohat(:,t+1))>1)
        [u,v]=max(crit);
        
        Ohat(:,t+1)=0;
        Ohat(v,t+1)=1;
    end
%     x(t+1)=x(t)+s(t)*dt;
    if t ~= 1
    shat(t) = (xhat(t)-xhat(t-1))/dt;
    end
    xhat(t+1)=(1-lambdaD*dt)*xhat(t)+Gamma*O(:,t);
    xhatOhat(t+1)=(1-lambdaD*dt)*xhatOhat(t)+Gamma*Ohat(:,t);
end

figure; hold on
for i = 1:N
line(xtime,-Ohat(i,:)+i,'color','k'); 
end
figure; hold on
for i = 1:N
line(xtime,-O(i,:)+i,'color','k'); 
end
% 
% figure;
% hold on
% plot(r')
% 
% figure;
% hold on
% plot(V')

figure;
hold on
% line(xtime(2:end),(diff(xhat)/dt)')
line(xtime,shat')

figure;
hold on
line(xtime,xhat./max(xhat),'color','k')
line(xtime,xhatOhat./max(xhatOhat),'color','r')
