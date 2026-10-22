function jequaloneDist=HuggettVenturaYaron2011_jequaloneDistFn(a_grid,z_grid,n_a,n_z,mean_logh1,mean_logability,stddev_logh1,stddev_logability,FTcorr_logh1logability,mu_z,sigma_z,ability_grid,N_i)
% Initial distribution: zero assets, and (h1,ability) bivariate log-normal, HVY2011 Table 2.
% Note: to get ability_grid and N_i we have to add them as inputs
%
% As in the authors' code (hcapshocks-ge.f), human capital at age 23 is exp(z1)*h1: the initial
% human capital h1 drawn from the log-normal is hit by a human capital shock z1~N(mu_z,sigma_z^2)
% before the first period. So log human capital at age 23 is normal with mean mean_logh1+mu_z and
% variance stddev_logh1^2+sigma_z^2, and the same covariance with log ability as log(h1).
% [(mean_logh1,stddev_logh1) are the parameters of h1, which is what HVY2011 Table 2 reports.]

% We parametrize the Covar matrix of (h1,a) using method of AH2021: https://github.com/robertdkirkby/ParametrizeCoVarianceMatrix
[CorrMatrix, ~ ] = GFT_inverse_mapping(FTcorr_logh1logability, 10^(-9));
corr_logh1logability=CorrMatrix(1,2);

asset_grid=a_grid(1:n_a(1));
h_grid=a_grid(n_a(1)+1:end);

logh1loga_Mean=[mean_logh1+mu_z, mean_logability];
logh1loga_CoVarMatrix=[stddev_logh1^2+sigma_z^2, corr_logh1logability*stddev_logh1*stddev_logability; corr_logh1logability*stddev_logh1*stddev_logability, stddev_logability^2];

% Now we can find the probabilities of a bivariate log-normal distribution over the existing grids
P=gpuArray(MVNormal_ProbabilitiesOnGrid([gather(log(h_grid));gather(log(ability_grid))],gather(logh1loga_Mean), gather(logh1loga_CoVarMatrix), [n_a(2),N_i])); % note: n_a(2) is number of points for human capital
% Note: P is on defined on (h,ability).

% NOTE: The need for all the gather() and gpuArray() in the line above is
% because there is an obscure error in mvncdf() on the gpu. It will likely
% be fixed (and so all the gather() and gpuArray() can be removed) in a
% future version of matlab (this was written late 2023).
% Explained at: https://au.mathworks.com/matlabcentral/answers/2056774-obscure-error-in-mvncdf

jequaloneDist=zeros([n_a,N_i],'gpuArray'); % Put no households anywhere on grid (I just leave this here to show what size it is)

% Initial distribution of agents at birth (j=1): zero assets
[~,zeroassetindex]=min(abs(asset_grid)); % the asset grid includes negative values (borrowing), and zero
jequaloneDist(zeroassetindex,:,:)=shiftdim(P,-1); % joint log-normal distribution onto our existing grids

end
