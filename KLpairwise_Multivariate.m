function kl_pairwise = KLpairwise_Multivariate(m , v)

% input: m (mean) as a row vector and v (var) as a [mxm] matrix estimate for each stim from KLprep_Multivariate

nstims = size(m,1);

kl_pairwise = nan(nstims,nstims);
for istim1 = 1:nstims
    for istim2 = istim1:nstims
    
        kl = KL_Multivar_symmetric(m(istim1,:)', squeeze(v(istim1,:,:)), m(istim2,:)', squeeze(v(istim2,:,:)));
        kl_pairwise(istim1,istim2) = kl;
        
    end
    
    fileID = fopen('matlog.txt','w');
    t = datestr(datetime('now'));
    fprintf(fileID,'%i %s\n',istim1,t);
    fclose(fileID);
end

idx = true(nstims);
idx = ~tril(idx);

kl_pairwise = kl_pairwise(idx);
