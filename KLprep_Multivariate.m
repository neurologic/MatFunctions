function [m, v, s] = KLprep_Multivariate(data)

% data is a [ stims x estimates (trials) x dimensions (time samples) ] mat

nstims = size(data,1);
nestimates = size(data,2);
ndimensions = size(data,3);

m = nan(nstims,ndimensions); % ends up as [stim x dimensions(samps)] size
v = nan(nstims,ndimensions,ndimensions) ;% ends up as [stim x dimensions(samps) x dimensions(samps)] size
for istim = 1:nstims
   thisstim = squeeze(data(istim,:,:));
   m(istim,:) = mean(thisstim,1);
   [sigma,shrinkage]=covCor(thisstim);
   v(istim,:,:) = sigma;
   s(istim) = shrinkage;
end
