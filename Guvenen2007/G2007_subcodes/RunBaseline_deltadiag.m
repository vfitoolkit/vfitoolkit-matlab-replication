% RunBaseline.m
% Stage (c) of the Guvenen (2007 AER) replication: the HIP life-cycle model with
% Bayesian learning about the income profile, for lambda in {0, 0.4, 0.62, 1}.
% Requires GPU (20GB: pi_z_J alone is ~8.8GB at the default grid). Run from baseline/.
%
% KNOWN-ALPHA DESIGN (see NOTES.md stage (c) for the full story): agents know their
% intercept alpha_i (the paper's footnote 4: alpha-uncertainty "turns out not to play
% an important role"; learning about it resolves within a few years). alpha is a
% 5-point discrete type; the filter runs on y - alpha = beta*exper + z + eps, so the
% belief chain (betahat, zhat, v) is TYPE-INVARIANT and is built once per lambda; the
% types differ only by the income-scale parameter alpha_i, a shift of the betahat grid
% (betabar_i = betabar + (sab/s2a)*(alpha_i-alphabar), preserving cov(alpha,beta)), and
% get solved as 5 separate runs aggregated by the law of total variance.
% Why not the 4-state chain of the brief: posterior means are martingales, and lottery
% discretization adds first-order variance per period; measured errors at feasible
% 4-dim grids were 65-133%. In 3-dim, nb=61 points on betahat gets ~1.6%.
%
% lambda semantics: knowing alpha reveals R2 = sab^2/(s2a*s2b) ~ 2.4% of beta-variance;
% the private prior knows fraction lambda of the REMAINING (conditional) variance:
% prior var = (1-lambda)*s2b_cond, s2b_cond = s2b*(1-R2).
%
% delta is recalibrated to W/Y=4 (total income incl r*a) for each lambda (paper p.700).
% Checks are inline, PASS/FAIL to the diary. A FAIL is a report, not a retune.

% ============================================================================
% DIAGNOSTIC COPY of RunBaseline.m (NOT the deliverable). Pins delta=0.966 (the
% paper's stated HIP value) for ALL lambda, SKIPS the W/Y=4 calibration loop, and
% reports BOTH the var(log c) rises AND the resulting W/Y (total-income and
% labour-only). Tests whether our accurate-grid model can match the paper's Fig-8
% rises and its W/Y=4 target at the same time, or whether delta=0.966 and W/Y=4 are
% mutually inconsistent once the grid is accurate (hypothesis: the paper's coarse
% 12-pt wealth grid under-measures wealth, so its calibration lands high at 0.966).
% Calibration loop + figures commented out for speed (~20 solves, ~5 min).
% ============================================================================
clearvars -except doPart
if exist('RunBaseline_deltadiag_diary.txt','file'); delete('RunBaseline_deltadiag_diary.txt'); end
diary('RunBaseline_deltadiag_diary.txt');
fprintf('=== RunBaseline_deltadiag.m (delta=0.966 fixed) run %s ===\n',char(datetime('now')));
npass=0; nfail=0;

%% Parameters
Jwork=40; Jret=31; N_j=Jwork+Jret; agevec=25:95;

% Income process, UNCONDITIONAL (Table 1 row 2, HIP; sigma_ab per NOTES.md)
s2a_u=0.022; s2b_u=0.00038; sab_u=-0.00045;
Params.rho=0.821;
Params.sigma2_eta=0.029;
Params.sigma2_eps=0.047;
Params.betabar=0.009;
Params.alphabar=1.5;
Params.exper1=0;
Params.T=Jwork;
% Filter (conditional on knowing alpha): zero alpha-uncertainty, conditional beta variance
R2ab=sab_u^2/(s2a_u*s2b_u);
Params.sigma2_alpha=0;
Params.sigma_ab=0;
Params.sigma2_beta=s2b_u*(1-R2ab);
fprintf('R2 of beta on alpha = %.4f; conditional s2b = %.6f (unconditional %.6f)\n',R2ab,Params.sigma2_beta,s2b_u);

% Preferences etc (Table 3)
Params.crra=2;
Params.Pb=0.96; Params.r=1/Params.Pb-1;
Params.delta=0.966;          % paper's HIP value; recalibrated per lambda below
Params.Phibar=1/1.40;

Params.exper_j=[0:Jwork-1,(Jwork-1)*ones(1,Jret)];
Params.exper_T=Jwork-1;
Params.workret_j=[ones(1,Jwork),zeros(1,Jret)];
Params.mewj=ones(1,N_j)/N_j;
AgeWeightParamNames={'mewj'};
DiscountFactorParamNames={'delta'};

lambdavec=[0,0.4,0.62,1];
ilam_base=3;
do_deltacalib=0;             % DIAGNOSTIC: skip calibration, use fixed Params.delta=0.966 (set above) for all lambda

%% Grid sizes
n_d=0; d_grid=[];
n_z=[61,13,5];               % [betahat, zhat, v]; chain errors at this size: betahat-var
                             % ~1.6%, var(log y) ~2.3% (measured; see NOTES.md).
                             % Refinement option [61,15,5] costs pi_z_J ~11.7GB on GPU.
N_z=prod(n_z); nb=n_z(1);
Params.nSigmasBelief=3;
discopts.parallel=0;
discopts.nSigmas=2.5;        % paper truncates the income distribution at 2.5 sd (p.702)
fprintf('n_z=[%d,%d,%d], N_z=%d; pi_z_J will need %.1f GB\n',n_z(1),n_z(2),n_z(3),N_z,N_z^2*(N_j-1)*8/1e9);

%% alpha types (known intercepts; moment-matched 5-point discretization)
n_types=5;
[alphagrid,w_types]=discretizeIIDNormal_TanakaToda(Params.alphabar,sqrt(s2a_u),n_types,discopts);
alphagrid=gather(alphagrid); w_types=gather(w_types);
dbet=(sab_u/s2a_u)*(alphagrid-Params.alphabar); % betabar tilt per type: preserves cov(alpha,beta)
va_disc=w_types'*(alphagrid.^2)-(w_types'*alphagrid)^2;
cab_disc=w_types'*(alphagrid.*(Params.betabar+dbet))-(w_types'*alphagrid)*(w_types'*(Params.betabar+dbet));
fprintf('alpha types: discrete var(alpha)=%.4f (target %.4f), cov(alpha,beta)=%.6f (target %.6f)\n',va_disc,s2a_u,cab_disc,sab_u);
if abs(va_disc/s2a_u-1)<1e-6 && abs(cab_disc/sab_u-1)<1e-6; npass=npass+1;
else; nfail=nfail+1; fprintf('FAIL: discrete alpha-type moments\n'); end

%% Natural borrowing limit (public prior = UNCONDITIONAL parameters, p.702; common to all lambda)
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
fprintf('Natural borrowing limit: min(Wbar)=%.3f, Wbar at last working age=%.3f\n',min(Params.Wbar_j),Params.Wbar_j(Jwork));

%% Ybar_T and the asset grid, from the BASELINE-lambda chain (pension rule fixed across lambda)
Params.lambda=lambdavec(ilam_base);
kf=KalmanSetup(Params);
[z_grid_J,pi_z_J,jequaloneDistz,grids]=BuildPiZ(kf,Params,n_z,N_j,discopts);
pj=jequaloneDistz';
for jj=1:Jwork-1; pj=pj*pi_z_J(:,:,jj); end
ynode_T=grids.bgrid_J(grids.IB,Jwork)*Params.exper_T+grids.zgrid_J(grids.IZ,Jwork)+grids.vgrid_J(grids.IV,Jwork);
E_common=pj*exp(ynode_T); % chain part, betahat centered on betabar
Params.Ybar_T=E_common*(w_types'*exp(alphagrid+dbet*Params.exper_T));
fprintf('Ybar_T (mean income at last working age) = %.4f\n',Params.Ybar_T);
amin=min(Params.Wbar_j);
amax=60*Params.Ybar_T;
a_grid=unique([linspace(amin,0,97)'; (amax*(linspace(0,1,305).^3))']);
n_a=length(a_grid);
fprintf('Asset grid: %d points on [%.2f, %.2f]\n',n_a,amin,amax);
clear pi_z_J % rebuilt per lambda below

%% FnsToEvaluate (alpha_i is a per-type parameter; log guard per NOTES.md stage (b))
FnsToEvaluate.cons=@(aprime,a,bh,zh,v,workret_j,exper_j,Pb,Phibar,Ybar_T,exper_T,alpha_i) ...
    a-Pb*aprime+workret_j*exp(alpha_i+bh*exper_j+zh+v)+(1-workret_j)*Phibar*Ybar_T*( ...
    (exp(alpha_i+bh*exper_T+zh+v)/Ybar_T<0.3)*(0.9*exp(alpha_i+bh*exper_T+zh+v)/Ybar_T) ...
    +(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T>=0.3)*(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T<=2)*(0.27+0.32*(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T-0.3)) ...
    +(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T>2)*(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T<=4.1)*(0.81+0.15*(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T-2)) ...
    +(exp(alpha_i+bh*exper_T+zh+v)/Ybar_T>4.1)*1.1);
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
FnsToEvaluate.logy=@(aprime,a,bh,zh,v,workret_j,exper_j,exper_T,alpha_i) ...
    workret_j*(alpha_i+bh*exper_j+zh+v)+(1-workret_j)*(alpha_i+bh*exper_T+zh+v);
FnsCalib.assets=FnsToEvaluate.assets;
FnsCalib.income=FnsToEvaluate.income;

ReturnFn=@(aprime,a,bh,zh,v,workret_j,exper_j,Wbar_j,Pb,crra,Phibar,Ybar_T,exper_T,alpha_i) ...
    ReturnFn_HIP(aprime,a,bh,zh,v,workret_j,exper_j,Wbar_j,Pb,crra,Phibar,Ybar_T,exper_T,alpha_i);
vfoptions.verbose=0;
vfoptions.divideandconquer=1; % ReturnFn matrix would be ~8GB done in one block
vfoptions.gridinterplayer=1;  % interpolate aprime between asset nodes (smooth policy)
vfoptions.ngridinterp=10;
simoptions.gridinterplayer=1; % MUST match vfoptions for StationaryDist/LifeCycleProfiles
simoptions.ngridinterp=10;


%% Loop over lambda
results=struct();
for ilam=1:length(lambdavec)
    Params.lambda=lambdavec(ilam);
    fprintf('\n============ lambda = %.2f ============\n',Params.lambda);
    kf=KalmanSetup(Params);
    [z_grid_J,pi_z_J,jequaloneDistz,grids]=BuildPiZ(kf,Params,n_z,N_j,discopts);

    %% Pre-solve chain checks
    maxrowdev=0;
    for jj=1:N_j-1
        maxrowdev=max(maxrowdev,max(abs(sum(pi_z_J(:,:,jj),2)-1)));
    end
    fprintf('max |row sum - 1| over pi_z_J slices = %.3e (tol 1e-10)\n',maxrowdev);
    if maxrowdev<1e-10; npass=npass+1; else; nfail=nfail+1; fprintf('FAIL: pi_z_J row sums\n'); end
    % chain moments vs Kalman theory (within type), plus AGGREGATE var(log y) vs the
    % unconditional Table-2 decomposition (adds the between-type part analytically)
    pj=jequaloneDistz';
    maxdrift=0; maxerrVb=0; maxerrVyagg=0;
    for t=1:Jwork
        if t>1; pj=pj*pi_z_J(:,:,t-1); end
        bnode=grids.bgrid_J(grids.IB,t);
        ynode=bnode*Params.exper_j(t)+grids.zgrid_J(grids.IZ,t)+grids.vgrid_J(grids.IV,t);
        Eb=pj*bnode; Vb=pj*(bnode.^2)-Eb^2;
        Ey=pj*ynode; Vy=pj*(ynode.^2)-Ey^2; % within-type
        Vb_th=Params.sigma2_beta-kf.Ppred(2,2,t);
        expr=Params.exper_j(t);
        Mtype=alphagrid+dbet*expr; % between-type mean log income (net of common part)
        Vbetween=w_types'*(Mtype.^2)-(w_types'*Mtype)^2;
        Vy_agg=Vy+Vbetween;
        Vy_th=s2a_u+Params.sigma2_eps+Params.sigma2_eta*(1-Params.rho^(2*t))/(1-Params.rho^2) ...
            +2*sab_u*expr+s2b_u*expr^2;
        maxdrift=max(maxdrift,abs(Eb-Params.betabar)/sqrt(Params.sigma2_beta));
        maxerrVb=max(maxerrVb,abs(Vb-Vb_th)/Params.sigma2_beta);
        maxerrVyagg=max(maxerrVyagg,abs(Vy_agg/Vy_th-1));
    end
    fprintf('chain vs theory: |E[betahat]-betabar|/sd = %.3e (tol 1e-2)\n',maxdrift);
    fprintf('chain vs theory: max |Var(betahat)-theory|/s2b = %.3e (tol 3e-2)\n',maxerrVb);
    fprintf('chain+types vs Table 2: max rel err var(log y) = %.3e (tol 3e-2)\n',maxerrVyagg);
    if maxdrift<1e-2; npass=npass+1; else; nfail=nfail+1; fprintf('FAIL: betahat drift\n'); end
    if maxerrVb<3e-2; npass=npass+1; else; nfail=nfail+1; fprintf('FAIL: Var(betahat) path (raise nb)\n'); end
    if maxerrVyagg<3e-2; npass=npass+1; else; nfail=nfail+1; fprintf('FAIL: aggregate var(log y) (raise nz/nb)\n'); end

    [~,a0ind]=min(abs(a_grid));
    jequaloneDist=zeros(n_a,N_z,'gpuArray');
    jequaloneDist(a0ind,:)=jequaloneDistz';
    jequaloneDist=reshape(jequaloneDist,[n_a,n_z]);

    %% delta calibration (aggregate W/Y over types = 5 solves per iteration)
    if do_deltacalib==1
        fprintf('--- Calibrating delta to W/Y=4 (lambda=%.2f) ---\n',Params.lambda);
        dlo=0.90; dhi=0.978;
        for calibit=1:8
            Params.delta=(dlo+dhi)/2;
            meana_agg=zeros(1,N_j); meany_agg=zeros(1,N_j);
            for ii=1:n_types
                Params.alpha_i=alphagrid(ii);
                zgJ_type=z_grid_J; zgJ_type(1:nb,:)=z_grid_J(1:nb,:)+dbet(ii);
                [~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,zgJ_type,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
                StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z_J,Params,simoptions);
                ACS=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsCalib,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,zgJ_type,simoptions);
                meana_agg=meana_agg+w_types(ii)*gather(ACS.assets.Mean);
                meany_agg=meany_agg+w_types(ii)*gather(ACS.income.Mean);
            end
            WYtot=sum(Params.mewj.*meana_agg)/sum(Params.mewj.*(meany_agg+Params.r*meana_agg));
            fprintf('  delta=%.5f: W/Y total-income=%.3f\n',Params.delta,WYtot);
            if WYtot>4; dhi=Params.delta; else; dlo=Params.delta; end
        end
        Params.delta=(dlo+dhi)/2;
        fprintf('calibrated delta(lambda=%.2f) = %.4f (paper reports 0.966 for HIP baseline)\n',Params.lambda,Params.delta);
    end

    %% Final solve per type, then aggregate by the law of total variance
    Mc=zeros(n_types,N_j); My=zeros(n_types,N_j); Ma=zeros(n_types,N_j);
    Mlc=zeros(n_types,N_j); Vlc=zeros(n_types,N_j); Mly=zeros(n_types,N_j); Vly=zeros(n_types,N_j);
    topmass=0; mono=true;
    for ii=1:n_types
        Params.alpha_i=alphagrid(ii);
        zgJ_type=z_grid_J; zgJ_type(1:nb,:)=z_grid_J(1:nb,:)+dbet(ii);
        tic;
        [V,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,zgJ_type,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
        vftime=toc;
        Vtest=reshape(gather(V),[n_a,N_z,N_j]);
        for jj=[1,Jwork,N_j]
            for zz=[1,ceil(N_z/2),N_z]
                vslice=Vtest(:,zz,jj); vslice=vslice(isfinite(vslice));
                mono=mono && all(diff(vslice)>=-1e-9);
            end
        end
        StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z_J,Params,simoptions);
        ACS=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,zgJ_type,simoptions);
        Mc(ii,:)=gather(ACS.cons.Mean); My(ii,:)=gather(ACS.income.Mean); Ma(ii,:)=gather(ACS.assets.Mean);
        Mlc(ii,:)=gather(ACS.logc.Mean); Vlc(ii,:)=gather(ACS.logc.Variance);
        Mly(ii,:)=gather(ACS.logy.Mean); Vly(ii,:)=gather(ACS.logy.Variance);
        topmass=topmass+w_types(ii)*gather(sum(StationaryDist(end-2:end,:,:,:,:),'all'));
        fprintf('type %d (alpha=%.3f): vf time %.1fs\n',ii,alphagrid(ii),vftime);
    end
    if mono; npass=npass+1; fprintf('V nondecreasing in assets (all types): PASS\n');
    else; nfail=nfail+1; fprintf('FAIL: V not monotone in assets\n'); end

    meanc=w_types'*Mc; meany=w_types'*My; meana=w_types'*Ma;
    varlogc=w_types'*(Vlc+Mlc.^2)-(w_types'*Mlc).^2; % law of total variance
    varlogy=w_types'*(Vly+Mly.^2)-(w_types'*Mly).^2;

    %% Per-lambda post-solve checks
    maxerr=0;
    for aa=[30 45 64]
        t=aa-24; expr=Params.exper_j(t);
        van=s2a_u+Params.sigma2_eps+Params.sigma2_eta*(1-Params.rho^(2*t))/(1-Params.rho^2)+2*sab_u*expr+s2b_u*expr^2;
        fprintf('var(log y) age %d: model %.4f, analytic %.4f\n',aa,varlogy(t),van);
        maxerr=max(maxerr,abs(varlogy(t)/van-1));
    end
    if maxerr<4e-2; npass=npass+1; else; nfail=nfail+1; fprintf('FAIL: var(log y) from dist vs analytic\n'); end
    WYtot=sum(Params.mewj.*meana)/sum(Params.mewj.*(meany+Params.r*meana));
    WYlab=sum(Params.mewj.*meana)/sum(Params.mewj.*meany); % labour+pension-only denominator
    fprintf('W/Y at delta=%.4f: total-income = %.2f, labour+pension = %.2f (target 4)\n',Params.delta,WYtot,WYlab);
    if do_deltacalib==0 || (WYtot>3.7 && WYtot<4.3); npass=npass+1; else; nfail=nfail+1; fprintf('FAIL: W/Y not at target\n'); end
    gret=mean(meanc(Jwork+5:N_j-1)./meanc(Jwork+4:N_j-2))-1;
    geuler=(Params.delta*(1+Params.r))^(1/Params.crra)-1;
    fprintf('retirement mean-c growth %.4f%%/yr vs Euler prediction %.4f%%/yr\n',100*gret,100*geuler);
    if abs(gret-geuler)<0.001; npass=npass+1; else; nfail=nfail+1; fprintf('FAIL: retirement Euler drift\n'); end
    fprintf('mass in top 3 asset points (type-weighted) = %.3e (tol 1e-5)\n',topmass);
    if topmass<1e-5; npass=npass+1; else; nfail=nfail+1; fprintf('FAIL: mass at top of asset grid\n'); end

    results(ilam).lambda=Params.lambda;
    results(ilam).delta=Params.delta;
    results(ilam).meanc=meanc; results(ilam).meany=meany; results(ilam).meana=meana;
    results(ilam).varlogc=varlogc; results(ilam).varlogy=varlogy;
    results(ilam).WYtot=WYtot;
    results(ilam).WYlab=WYlab;
    results(ilam).risevarlogc=varlogc(Jwork+1)-varlogc(1);
    fprintf('var(log c): %.3f at 25, %.3f at 65, RISE = %.3f\n',varlogc(1),varlogc(Jwork+1),results(ilam).risevarlogc);
end

%% Cross-lambda checks (the stage-(c) headline: Fig 8 / brief 1.8 check 7)
rises=[results.risevarlogc];
fprintf('\n--- Cross-lambda results ---\n');
fprintf('lambda:            ');fprintf(' %8.2f',lambdavec);fprintf('\n');
fprintf('delta*:            ');fprintf(' %8.4f',[results.delta]);fprintf('\n');
fprintf('var(log c) rise:   ');fprintf(' %8.3f',rises);fprintf('\n');
fprintf('W/Y total-income:  ');fprintf(' %8.2f',[results.WYtot]);fprintf('\n');
fprintf('W/Y labour+pens:   ');fprintf(' %8.2f',[results.WYlab]);fprintf('\n');
fprintf('(paper Fig 8 rises: lambda=0 ~0.32, lambda*=0.62 ~0.21 [matches US data], lambda=1 ~0.14; paper W/Y target = 4)\n');
fprintf('DIAGNOSTIC READING: if rises now ~=paper AND W/Y>>4, then delta=0.966 and W/Y=4 are inconsistent in our\n');
fprintf('  accurate-grid model (paper coarse grid under-measures wealth). If rises still high, delta is not the cause.\n');
if all(diff(rises)<0); fprintf('rise in var(log c) decreasing in lambda: PASS\n');
else; fprintf('NOTE: rise not monotone decreasing in lambda\n'); end

%% Figures (commented out for the fast diagnostic)
% fig=figure(1); clf; hold on;
% for ilam=1:length(lambdavec)
%     plot(25:65,results(ilam).varlogc(1:Jwork+1),'LineWidth',1.5);
% end
% xlabel('age'); ylabel('var(log c)'); legend('\lambda=0','\lambda=0.4','\lambda=0.62','\lambda=1','Location','northwest');
% title('Age-inequality profile of consumption by \lambda (cf. Fig 8)');
% saveas(fig,'SavedOutput/Graphs/RunBaseline_fig1_varlogc_bylambda.png');

%% Save and summarize
save('RunBaseline_deltadiag_results.mat','results','lambdavec','Params','a_grid','n_z','n_a','N_j','alphagrid','w_types','dbet','-v7.3');
fprintf('\n=== RunBaseline summary: %d passed, %d failed ===\n',npass,nfail);
diary off;
