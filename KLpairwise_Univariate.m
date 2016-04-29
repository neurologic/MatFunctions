function kl_pairwise = KLpairwise_Univariate(m , v)

% input: m (mean) and v (var) estimate for each stim from KLprep_Univariate

nstims = size(m,1);

kl_pairwise = nan(nstims,nstims);
for istim1 = 1:nstims
    for istim2 = istim1:nstims
    
        kl = KL_Univar_symmetric(m(istim1), v(istim1), m(istim2), v(istim2));
        kl_pairwise(istim1,istim2) = kl;
        
    end
end

idx = true(nstims);
idx = ~tril(idx);

kl_pairwise = kl_pairwise(idx);
