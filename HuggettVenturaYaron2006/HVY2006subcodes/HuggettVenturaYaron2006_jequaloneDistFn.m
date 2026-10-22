function jequaloneDist=HuggettVenturaYaron2006_jequaloneDistFn(h_grid,z_grid,n_a,n_z,mean_logh1,mean_logability,stddev_logh1,stddev_logability,FTcorr_logh1logability, ability_grid,N_i)
% Note: to get ability_grid and N_i we have to add them as inputs

% We parametrize the Covar matrix of (h1,a) using method of AH2021: https://github.com/robertdkirkby/ParametrizeCoVarianceMatrix
[CorrMatrix, ~ ] = GFT_inverse_mapping(FTcorr_logh1logability, 10^(-9));
logh1loga_CoVarMatrix = corr2cov_homemade([stddev_logh1,stddev_logability],CorrMatrix);


logh1loga_Mean=[mean_logh1,mean_logability];

% Now we can find the probabilities of a bivariate log-normal distribution over the existing grids
P=gpuArray(MVNormal_ProbabilitiesOnGrid([gather(log(h_grid));gather(log(ability_grid))],gather(logh1loga_Mean), gather(logh1loga_CoVarMatrix), [n_a,N_i])); % note: n_a is number of points for human capital
% Note: P is on defined on (h,ability).

% NOTE: The need for all the gather() and gpuArray() in the line above is
% because there is an obscure error in mvncdf() on the gpu. It will likely
% be fixed (and so all the gather() and gpuArray() can be removed) in a 
% future version of matlab (this was written late 2023).
% Explained at: https://au.mathworks.com/matlabcentral/answers/2056774-obscure-error-in-mvncdf

% Initial distribution of agents at birth (j=1)
jequaloneDist=P; % joint log-normal distribution onto our existing grids



end