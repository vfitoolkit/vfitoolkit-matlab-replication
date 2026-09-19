% ============================================================================
% RunBaseline_gridtest.m  --  DIAGNOSTIC (not the deliverable).
%
% Question: is the paper's reported delta=0.966-at-W/Y=4 a WEALTH-GRID ARTIFACT?
% Guvenen (2007) solves the value function on a 12-point wealth grid closed by a
% global polynomial regression; we use 501 points with direct (divide-and-conquer)
% VFI. At the paper's stated delta=0.966 our fine grid gives W/Y ~= 5.7, not 4.
%
% Test: hold delta=0.966 and EVERYTHING else fixed, and sweep ONLY the number of
% wealth-grid points from 501 down to Guvenen's 12. If measured W/Y falls toward 4
% as the grid coarsens, then a coarse grid under-measures wealth -- i.e. the paper,
% calibrating on a coarse grid, would have seen W/Y=4 at delta=0.966 while the
% accurate answer at that delta is ~5.7. That both explains the delta gap
% (0.966 coarse vs ~0.95 accurate, for the same W/Y=4 target) and confirms the
% accurate solution is the one to report.
%
% Caveat: our coarse grid uses piecewise (DC) interpolation, whereas Guvenen closes
% his 12 points with a smooth global polynomial, so this reproduces the DIRECTION of
% the grid effect, not his exact pipeline. Baseline lambda=0.62 only, for speed.
% ============================================================================
clearvars -except doPart
if exist('RunBaseline_gridtest_diary.txt','file'); delete('RunBaseline_gridtest_diary.txt'); end
diary('RunBaseline_gridtest_diary.txt');
fprintf('=== RunBaseline_gridtest.m (delta=0.966, wealth-grid sweep, lambda=0.62) run %s ===\n',char(datetime('now')));

%% Parameters (identical to RunBaseline; delta fixed at the paper's 0.966)
Jwork=40; Jret=31; N_j=Jwork+Jret; agevec=25:95;
s2a_u=0.022; s2b_u=0.00038; sab_u=-0.00045;
Params.rho=0.821; Params.sigma2_eta=0.029; Params.sigma2_eps=0.047;
Params.betabar=0.009; Params.alphabar=1.5; Params.exper1=0; Params.T=Jwork;
R2ab=sab_u^2/(s2a_u*s2b_u);
Params.sigma2_alpha=0; Params.sigma_ab=0; Params.sigma2_beta=s2b_u*(1-R2ab);
Params.crra=2; Params.Pb=0.96; Params.r=1/Params.Pb-1;
Params.delta=0.966;          % the paper's stated HIP value, held fixed for this test
Params.Phibar=1/1.40;
Params.exper_j=[0:Jwork-1,(Jwork-1)*ones(1,Jret)];
Params.exper_T=Jwork-1;
Params.workret_j=[ones(1,Jwork),zeros(1,Jret)];
Params.mewj=ones(1,N_j)/N_j;
AgeWeightParamNames={'mewj'};
DiscountFactorParamNames={'delta'};
lambda=0.62; Params.lambda=lambda;

n_d=0; d_grid=[];
n_z=[61,13,5]; N_z=prod(n_z); nb=n_z(1);
Params.nSigmasBelief=3;
discopts.parallel=0; discopts.nSigmas=2.5;

% wealth-grid sizes to sweep (Guvenen used 12; ours is 501)
ngrid_list=[12, 25, 50, 100, 200, 301]; % top capped for GI memory at N_z=3965

%% alpha types (same as RunBaseline)
n_types=5;
[alphagrid,w_types]=discretizeIIDNormal_TanakaToda(Params.alphabar,sqrt(s2a_u),n_types,discopts);
alphagrid=gather(alphagrid); w_types=gather(w_types);
dbet=(sab_u/s2a_u)*(alphagrid-Params.alphabar);

%% Natural borrowing limit (same as RunBaseline)
ymin=zeros(1,N_j); zmin=0;
for jj=1:Jwork
    zmin=Params.rho*zmin-2.5*sqrt(Params.sigma2_eta);
    ymin(jj)=exp((Params.alphabar-2.5*sqrt(s2a_u))+(Params.betabar-2.5*sqrt(s2b_u))*Params.exper_j(jj)+zmin-2.5*sqrt(Params.sigma2_eps));
end
ymin(Jwork+1:N_j)=Params.Phibar*0.9*ymin(Jwork);
M=zeros(1,N_j+1);
for jj=N_j:-1:1
    M(jj)=ymin(jj)+Params.Pb*M(jj+1);
end
Params.Wbar_j=[-M(2:N_j),0];

%% Belief chain (built once; only the wealth grid changes below)
kf=KalmanSetup(Params);
[z_grid_J,pi_z_J,jequaloneDistz,grids]=BuildPiZ(kf,Params,n_z,N_j,discopts);
pj=jequaloneDistz';
for jj=1:Jwork-1; pj=pj*pi_z_J(:,:,jj); end
ynode_T=grids.bgrid_J(grids.IB,Jwork)*Params.exper_T+grids.zgrid_J(grids.IZ,Jwork)+grids.vgrid_J(grids.IV,Jwork);
Params.Ybar_T=(pj*exp(ynode_T))*(w_types'*exp(alphagrid+dbet*Params.exper_T));
amin=min(Params.Wbar_j);
amax=60*Params.Ybar_T;
fprintf('lambda=%.2f, delta=%.4f, Ybar_T=%.3f, asset range [%.2f, %.2f]\n',lambda,Params.delta,Params.Ybar_T,amin,amax);

%% FnsToEvaluate (same as RunBaseline; only assets+income needed here)
FnsToEvaluate.income=@(aprime,a,bh,zh,v,workret_j,exper_j,Phibar,Ybar_T,exper_T,alpha_i) ...
    workret_j*exp(alpha_i+bh*exper_j+zh+v)+(1-workret_j)*Phibar*Ybar_T*( ...
    (exp(alpha_i+bh*exper_T+zh+v)/Ybar_T<0.3)*(0.9*exp(alpha_i+bh*exper_T+zh+v)/Ybar_T) ...
    +(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T>=0.3)*(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T<=2)*(0.27+0.32*(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T-0.3)) ...
    +(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T>2)*(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T<=4.1)*(0.81+0.15*(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T-2)) ...
    +(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T>4.1)*1.1);
FnsToEvaluate.assets=@(aprime,a,bh,zh,v) a;
FnsToEvaluate.logc=@(aprime,a,bh,zh,v,workret_j,exper_j,Pb,Phibar,Ybar_T,exper_T,alpha_i) ...
    log(max(1e-12,a-Pb*aprime+workret_j*exp(alpha_i+bh*exper_j+zh+v)+(1-workret_j)*Phibar*Ybar_T*( ...
    (exp(alpha_i+bh*exper_T+zh+v)/Ybar_T<0.3)*(0.9*exp(alpha_i+bh*exper_T+zh+v)/Ybar_T) ...
    +(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T>=0.3)*(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T<=2)*(0.27+0.32*(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T-0.3)) ...
    +(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T>2)*(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T<=4.1)*(0.81+0.15*(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T-2)) ...
    +(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T>4.1)*1.1)));

ReturnFn=@(aprime,a,bh,zh,v,workret_j,exper_j,Wbar_j,Pb,crra,Phibar,Ybar_T,exper_T,alpha_i) ...
    ReturnFn_HIP(aprime,a,bh,zh,v,workret_j,exper_j,Wbar_j,Pb,crra,Phibar,Ybar_T,exper_T,alpha_i);
vfoptions.verbose=0;
vfoptions.divideandconquer=1;
vfoptions.gridinterplayer=1;  % interpolate aprime between asset nodes (smooth policy)
vfoptions.ngridinterp=10;
simoptions.gridinterplayer=1; % MUST match vfoptions for StationaryDist/LifeCycleProfiles
simoptions.ngridinterp=10;

%% Sweep wealth-grid resolution
ng=length(ngrid_list);
res_na=zeros(1,ng); res_WYtot=zeros(1,ng); res_WYlab=zeros(1,ng);
res_rise=zeros(1,ng); res_top=zeros(1,ng); res_vft=zeros(1,ng);
for ig=1:ng
    ntot=ngrid_list(ig);
    % same two-piece construction as RunBaseline (even negative + cubic positive),
    % split ~24%/76% and scaled to ntot; 0 shared between the two pieces
    nneg=max(3,round(0.24*(ntot+1)));
    npos=(ntot+1)-nneg;
    a_grid=unique([linspace(amin,0,nneg)'; (amax*(linspace(0,1,npos).^3))']);
    n_a=length(a_grid);
    [~,a0ind]=min(abs(a_grid));
    jequaloneDist=zeros(n_a,N_z,'gpuArray');
    jequaloneDist(a0ind,:)=jequaloneDistz';
    jequaloneDist=reshape(jequaloneDist,[n_a,n_z]);

    Ma=zeros(n_types,N_j); My=zeros(n_types,N_j);
    Mlc=zeros(n_types,N_j); Vlc=zeros(n_types,N_j);
    topmass=0; tic;
    for ii=1:n_types
        Params.alpha_i=alphagrid(ii);
        zgJ_type=z_grid_J; zgJ_type(1:nb,:)=z_grid_J(1:nb,:)+dbet(ii);
        [~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,zgJ_type,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
        StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z_J,Params,simoptions);
        ACS=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,zgJ_type,simoptions);
        Ma(ii,:)=gather(ACS.assets.Mean); My(ii,:)=gather(ACS.income.Mean);
        Mlc(ii,:)=gather(ACS.logc.Mean); Vlc(ii,:)=gather(ACS.logc.Variance);
        topmass=topmass+w_types(ii)*gather(sum(StationaryDist(end-2:end,:,:,:,:),'all'));
    end
    vft=toc;
    meana=w_types'*Ma; meany=w_types'*My;
    varlogc=w_types'*(Vlc+Mlc.^2)-(w_types'*Mlc).^2;
    WYtot=sum(Params.mewj.*meana)/sum(Params.mewj.*(meany+Params.r*meana));
    WYlab=sum(Params.mewj.*meana)/sum(Params.mewj.*meany);
    rise=varlogc(Jwork+1)-varlogc(1);
    res_na(ig)=n_a; res_WYtot(ig)=WYtot; res_WYlab(ig)=WYlab;
    res_rise(ig)=rise; res_top(ig)=topmass; res_vft(ig)=vft;
    fprintf('n_a=%3d (target %3d): W/Y total=%.2f, labour=%.2f, var(logc) rise=%.3f, top-mass=%.2e, %.1fs\n',...
        n_a,ntot,WYtot,WYlab,rise,topmass,vft);
end

%% Summary
fprintf('\n--- Wealth-grid resolution sweep (lambda=0.62, delta=0.966 fixed) ---\n');
fprintf('n_a (Guvenen=12):   ');fprintf(' %8d',res_na);fprintf('\n');
fprintf('W/Y total-income:   ');fprintf(' %8.2f',res_WYtot);fprintf('\n');
fprintf('W/Y labour+pension: ');fprintf(' %8.2f',res_WYlab);fprintf('\n');
fprintf('var(log c) rise:    ');fprintf(' %8.3f',res_rise);fprintf('\n');
fprintf('(paper target W/Y=4; accurate 501-pt gives ~5.7 total at delta=0.966)\n');
fprintf('READING: if W/Y total falls toward 4 as n_a->12, a coarse grid under-measures\n');
fprintf('  wealth -> the paper''s W/Y=4-at-delta=0.966 is a grid artifact; ~0.95 (our value) is accurate.\n');

save('RunBaseline_gridtest_results.mat','ngrid_list','res_na','res_WYtot','res_WYlab','res_rise','res_top','res_vft','Params');
fprintf('\n=== RunBaseline_gridtest done ===\n');
diary off;
