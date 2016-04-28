function [m, v] = KLprep_Univariate(data)

% data is a [ stims x estimates (trials) x dimensions (time samples) ] mat

nstims = size(data,1);
nestimates = size(data,2);
ndimensions = size(data,3);

for istim = 1:nstims
   thisstim = squeeze(data(istim,:,:));
   thisstim = mean(thisstim,2);
   m(istim,1) = mean(thisstim);
   v(istim,1) = var(thisstim);
end