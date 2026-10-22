function jequaloneDist=HuggettVenturaYaron2011_jequaloneDistFn_wealth(a_grid,z_grid,n_a,n_z,mean_logh1,mean_logability,stddev_logh1,stddev_logability,FTw_ha,FTw_hk,FTw_ak,mean_logwealth,stddev_logwealth,wealthshift,earn23,mu_z,sigma_z,ability_grid,N_i)
% Initial distribution of the model with initial wealth differences (HVY2011 Section V.A): (log h1, log a, x) is trivariate
% normal, and initial assets are k1=(exp(x)-wealthshift)*earn23, a shifted log-normal in units of mean earnings at age 23
% (as in the authors' code distrib-estimation-initw.f: x~N(-0.277,0.921^2), shift 0.381, times data mean earnings at age 23).
% As in the benchmark, human capital at age 23 is exp(z1)*h1 (see HuggettVenturaYaron2011_jequaloneDistFn).
% The correlation matrix of (log h1, log a, x) is parametrized as in Archakov & Hansen (2021), by (FTw_ha, FTw_hk, FTw_ak).
%
% Rather than a trivariate grid, use that x conditional on (log h, log a) is normal: the probabilities of (h, a) are on the grids
% (as in the benchmark), and those of x given (h,a) are on nodes for x; each node gives a k1 that is put on the asset grid by
% linear interpolation between the two nearest grid points.

[CorrMatrix, ~ ] = GFT_inverse_mapping([FTw_ha;FTw_hk;FTw_ak], 10^(-9)); % lower triangle in column order: (h,a), (h,x), (a,x)
stddevvec=[stddev_logh1; stddev_logability; stddev_logwealth];
CoVarMatrix=CorrMatrix.*(stddevvec*stddevvec');

asset_grid=gather(a_grid(1:n_a(1)));
h_grid=gather(a_grid(n_a(1)+1:end));
ability_grid=gather(ability_grid);

% (log h at age 23, log a): h at age 23 is exp(z1)*h1, with z1 independent of everything else
logha_Mean=[mean_logh1+mu_z, mean_logability];
logha_CoVarMatrix=CoVarMatrix(1:2,1:2);
logha_CoVarMatrix(1,1)=logha_CoVarMatrix(1,1)+sigma_z^2;
P=gather(MVNormal_ProbabilitiesOnGrid([log(h_grid);log(ability_grid)],logha_Mean, logha_CoVarMatrix, [n_a(2),N_i])); % [n_a(2),N_i]

% x given (log h at age 23, log a): normal, mean linear in (log h, log a), constant variance
% (the covariance of x with log h at age 23 is that with log h1, as z1 is independent)
covx=CoVarMatrix(3,1:2);
Bx=covx/logha_CoVarMatrix;
condsd=sqrt(max(CoVarMatrix(3,3)-Bx*covx',10^(-12)));
[LogH,LogA]=ndgrid(log(h_grid),log(ability_grid));
condmean=mean_logwealth+Bx(1)*(LogH-logha_Mean(1))+Bx(2)*(LogA-logha_Mean(2)); % [n_a(2),N_i]

% Nodes for x, and the asset grid points (and weights) for the k1 of each node
nx=21;
xnodes=mean_logwealth+stddev_logwealth*linspace(-3,3,nx)';
xedges=[-Inf; (xnodes(1:end-1)+xnodes(2:end))/2; Inf];
k1=(exp(xnodes)-wealthshift)*earn23;
lowerindex=discretize(k1,asset_grid);
weightupper=(k1-asset_grid(lowerindex))./(asset_grid(lowerindex+1)-asset_grid(lowerindex));

jequaloneDist=zeros([n_a,N_i]);
for kk=1:nx
    probx=0.5*erfc(-(xedges(kk+1)-condmean)/(condsd*sqrt(2)))-0.5*erfc(-(xedges(kk)-condmean)/(condsd*sqrt(2))); % [n_a(2),N_i]
    mass=reshape(P.*probx,[1,n_a(2),N_i]);
    jequaloneDist(lowerindex(kk),:,:)=jequaloneDist(lowerindex(kk),:,:)+(1-weightupper(kk))*mass;
    jequaloneDist(lowerindex(kk)+1,:,:)=jequaloneDist(lowerindex(kk)+1,:,:)+weightupper(kk)*mass;
end
jequaloneDist=gpuArray(jequaloneDist/sum(jequaloneDist(:)));

end
