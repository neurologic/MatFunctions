function cutdata = removeSpikes(data,spikethresh,dt)
highpassdata = HighpassGeneral(data,1/dt);
highpassdata = highpassdata(:,2:end-2); %end before end so that don't count last transient as spike
% [spikesmat, gausstosmooth, spiketimes]=getspikesmat(highpassdata,spikethresh,dt);
spikesmat=[];
a=highpassdata>=spikethresh;
adif=diff(a')';
spikesmat=adif==1;
% [spikesmat]=getspikesmat(highpassdata,spikethresh,dt);
spiketimes = [];
for itrial=1:size(spikesmat,1)
    spiketimes{itrial}=find(spikesmat(itrial,:));%*dt;
end
cutdata = [];
for itrial=1:size(spikesmat,1)
    spks_trial = spiketimes{itrial};
    thistrial = data(itrial,:);
    if ~isempty(spks_trial)
        for ispike = 1:size(spks_trial,2)
            t1 = spks_trial(ispike)+5;
            
            spkwin_size = round((20/1000/dt)/2);
            if ispike == 1
                altbegin = [(t1 - spkwin_size),1];
            else altbegin = [(t1 - spkwin_size),1,spks_trial(ispike-1)+10];
            end
            
            if ispike == size(spks_trial,2);
                altend = [(t1 + spkwin_size),size(highpassdata,2),];
            else altend = [(t1 + spkwin_size),size(highpassdata,2),spks_trial(ispike+1)];
            end
            
            spk_win = [max(altbegin),min(altend)];
            
            tmpshape = thistrial(spk_win(1):spk_win(2));
            spkt = t1-spk_win(1);
            if (spkt-8) > 0 && (spkt+9) < size(tmpshape,2)
                peakt = min(find(tmpshape==max(tmpshape(spkt-8:spkt+9)))); %changed from spkt+20 to 10
            elseif (spkt-8) < 0 && (spkt+9) < size(tmpshape,2)
                peakt = min(find(tmpshape==max(tmpshape(1:spkt+9))));
            elseif (spkt-8) > 0 && (spkt+9) > size(tmpshape,2)
                peakt = min(find(tmpshape==max(tmpshape(spkt-8:end))));
            elseif (spkt-8) < 0 && (spkt+9) > size(tmpshape,2)
                peakt = min(find(tmpshape==max(tmpshape(1:end))));
            else
                return
            end
            dvdt = diff(tmpshape(1:peakt));
            if size(dvdt,2)>10
                dvdt_max = max(dvdt(1,end-10:end));
            else
                dvdt_max = max(dvdt);
            end
            peakind = max(find(dvdt == dvdt_max));
            %         dvdt_thr = 0.033 * dvdt_max;
            dvdt_thr = 0.1 * dvdt_max;
            
            spk_init = max(find(dvdt(1:peakind)<dvdt_thr));
            
            if isempty(spk_init)
                continue % BUT if this is catching actual spikes here that are part of bursts
                % i might want to instead cut these at sample 1
                % even if dvdt not below dvdt_max
            end
            
            spk_end = min(find(tmpshape(peakt:end)<tmpshape(spk_init)))+peakt-1;
            %if no value less than spk_init in that window)
            if isempty(spk_end)
                spk_end = min(find(tmpshape(peakt:end) == min(tmpshape(peakt:end)))) + peakt-1;
            end
            
            x = [spk_init,spk_end];
            xi = [spk_init:1:spk_end];
            y = [tmpshape(spk_init),tmpshape(spk_end)];
            yi = interp1(x,y,xi);
            
            cutshape = [tmpshape(1:spk_init-1),yi,tmpshape((spk_end+1):end)];
            
            startcut = spk_init + spk_win(1) -2;
            stopcut = spk_end + spk_win(1);
            thistrial = [thistrial(1:startcut),yi,data(itrial,stopcut:end)];
            cutdata(itrial,:) = thistrial;
        end
    end
    cutdata(itrial,:) = thistrial;
end