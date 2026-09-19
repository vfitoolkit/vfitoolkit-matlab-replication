function [z_grid_J,pi_z_J,jequaloneDistz,grids]=BuildPiZ(kf,Params,n_z,N_j,discopts)
% Build the exogenous belief-chain objects for the Guvenen (2007) HIP learning model,
% KNOWN-ALPHA version (see NOTES.md stage (c): alpha-learning dropped per the paper's
% footnote 4; alpha is a discrete type handled outside this function, and the filter
% runs on the alpha-free observation y - alpha = beta*exper + z + eps).
%
% Exogenous state: (betahat, zhat, v) with betahat=betahat_{t|t-1}, zhat=zhat_{t|t-1}
% the predictive beliefs and v_t ~ iid N(0,kf.Fv(t)) the income innovation (objectively
% iid by correct-prior Bayes). Within-type log income (net of alpha) is
%   ytilde_t = betahat*exper_t + zhat + v.
% Transition from working age t to t+1 (kf from KalmanSetup with sigma2_alpha=0,
% sigma_ab=0, sigma2_beta = the CONDITIONAL-on-alpha variance; kf.T = Jwork):
%   betahat' = betahat + K_beta,t*v     [martingale]
%   zhat'    = rho*(zhat + K_z,t*v)
%   v' ~ N(0,Fv_{t+1}) iid
% so pi_t = [2-dim bilinear lottery onto the age-(t+1) belief grid] tensor [phi_{t+1}(v')].
% From the last working period on, all transitions are identity (freeze).
%
% Grid sizes matter here (see NOTES.md): betahat is a slow martingale, and the lottery
% adds cross-sectional variance ~ increment*spacing per period. Measured: nb=61 keeps
% the betahat-variance error to ~1.6% (nb=41: ~5%; nb=7: ~120%). zhat is mean-reverting
% so nz~13 suffices. The chain-moment checks in RunBaseline verify per run.
%
% The chain is TYPE-INVARIANT: the transition and jequaloneDistz do not depend on the
% type's (alpha_i, betabar_i); RunBaseline shifts the betahat block of z_grid_J by
% (betabar_i - betabar) per type.
%
% Inputs: kf (KalmanSetup output), Params (betabar, lambda, rho, sigma2_beta
%   [conditional], sigma2_eta, nSigmasBelief), n_z=[nb,nz,nv], N_j, discopts
%   (for discretizeIIDNormal_TanakaToda; parallel=0, nSigmas=2.5 per the paper).
% Outputs: z_grid_J [sum(n_z),N_j] stacked; pi_z_J [prod(n_z),prod(n_z),N_j-1]
%   (NB: ~prod(n_z)^2*(N_j-1)*8 bytes -- 8.8GB at [61,13,5]); jequaloneDistz
%   [prod(n_z),1]; grids struct (per-age grids, weights, index expansions, sd paths).

Jwork=kf.T;
nb=n_z(1); nz=n_z(2); nv=n_z(3);
Nb=nb*nz; N_z=Nb*nv;
nSigB=Params.nSigmasBelief;
if mod(nb,2)~=1 || mod(nz,2)~=1
    error('BuildPiZ: nb and nz must be odd (period-1 point masses sit on the center grid point)')
end

%% Theory sd paths of the cross-section of predictive beliefs (within type)
sdb=zeros(1,Jwork); sdz=zeros(1,Jwork);
for t=1:Jwork
    sdb(t)=sqrt(max(Params.sigma2_beta-kf.Ppred(2,2,t),0));
    varz_t=Params.sigma2_eta*(1-Params.rho^(2*t))/(1-Params.rho^2); % population var of z_t (z_0=0)
    sdz(t)=sqrt(max(varz_t-kf.Ppred(3,3,t),0));
end
% floors keep period-1 (near-)degenerate grids strictly increasing for the lottery
sdb=max(sdb,1e-4*sqrt(Params.sigma2_beta));
sdz=max(sdz,1e-4*sqrt(Params.sigma2_eta));

%% Grids per working age
bgrid_J=Params.betabar+linspace(-nSigB,nSigB,nb)'.*sdb; % nb-by-Jwork
zgrid_J=               linspace(-nSigB,nSigB,nz)'.*sdz;
vgrid_J=zeros(nv,Jwork); wv_J=zeros(nv,Jwork);
for t=1:Jwork
    [vg,wv]=discretizeIIDNormal_TanakaToda(0,sqrt(kf.Fv(t)),nv,discopts);
    vgrid_J(:,t)=gather(vg); wv_J(:,t)=gather(wv);
end

% Period-1 betahat distribution: betahat_{1|0} = betabar + k, k ~ N(0,lambda*s2b_cond).
% For lambda>0 use a moment-matched TanakaToda discretization and adopt its grid as the
% period-1 betahat grid (initial variance exact); lambda=0 is a point mass.
if Params.lambda>0
    [bg1,wb1]=discretizeIIDNormal_TanakaToda(Params.betabar,sqrt(Params.lambda*Params.sigma2_beta),nb,discopts);
    bgrid_J(:,1)=gather(bg1); wb1=gather(wb1);
else
    wb1=zeros(nb,1); wb1((nb+1)/2)=1;
end

%% Stacked age-dependent grid; retirement reuses the last working-age grids (frozen)
z_grid_J=zeros(sum(n_z),N_j);
for jj=1:N_j
    t=min(jj,Jwork);
    z_grid_J(:,jj)=[bgrid_J(:,t);zgrid_J(:,t);vgrid_J(:,t)];
end

%% Transitions
% Joint ordering: first variable fastest (toolkit convention), linear index
% = ib + nb*(iz-1) + Nb*(iv-1); ndgrid is also first-dim-fastest.
[IB,IZ,IV]=ndgrid(1:nb,1:nz,1:nv);
IB=IB(:); IZ=IZ(:); IV=IV(:);
R=(1:N_z)';
pi_z_J=zeros(N_z,N_z,N_j-1);
for t=1:Jwork-1
    bp=bgrid_J(IB,t)+kf.K(2,t)*vgrid_J(IV,t);
    zp=Params.rho*(zgrid_J(IZ,t)+kf.K(3,t)*vgrid_J(IV,t));
    % 2-point bracket + weight per dimension of the age-(t+1) belief grid (clamp at edges)
    g=bgrid_J(:,t+1);
    x=min(max(bp,g(1)),g(nb)); ib1=min(discretize(x,g),nb-1); wb=(g(ib1+1)-x)./(g(ib1+1)-g(ib1));
    g=zgrid_J(:,t+1);
    x=min(max(zp,g(1)),g(nz)); iz1=min(discretize(x,g),nz-1); wz=(g(iz1+1)-x)./(g(iz1+1)-g(iz1));
    % 4 corners of the bilinear lottery
    cols=zeros(N_z,4); vals=zeros(N_z,4);
    k=0;
    for db=0:1
        for dz=0:1
            k=k+1;
            cols(:,k)=(ib1+db)+nb*((iz1+dz)-1);
            vals(:,k)=(db+(1-2*db)*wb).*(dz+(1-2*dz)*wz);
            % db=0 -> weight wb on the lower point; db=1 -> weight 1-wb on the upper
        end
    end
    B=sparse(repmat(R,4,1),cols(:),vals(:),N_z,Nb); % (belief,v) -> belief' lottery
    pi_z_J(:,:,t)=kron(wv_J(:,t+1)',full(B));       % tensor with next-age innovation dist
end
for jj=Jwork:N_j-1 % last-working-to-retirement and within retirement: freeze
    pi_z_J(:,:,jj)=eye(N_z);
end

%% Period-1 distribution over the joint z space
wz1=zeros(nz,1); wz1((nz+1)/2)=1; % zhat_1|0 = 0 (center)
jequaloneDistz=kron(wv_J(:,1),kron(wz1,wb1));

%% Pack grids for checks
grids.bgrid_J=bgrid_J; grids.zgrid_J=zgrid_J;
grids.vgrid_J=vgrid_J; grids.wv_J=wv_J;
grids.sdb=sdb; grids.sdz=sdz;
grids.IB=IB; grids.IZ=IZ; grids.IV=IV;

end
