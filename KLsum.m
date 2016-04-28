function kl_sums = KLsum(kl_pairwise,nsubsamps,nreps)


% nsubsamps = 1000;
% nreps = 1000;

allstimsize = size(kl_pairwise,1);

kl_sums = nan(1,nreps);

for irep = 1:nreps
   randrep = randperm(allstimsize);
   thesestims = randrep(1:nsubsamps);
   
   kl_sums(irep) = sum(kl_pairwise(thesestims));
   
end