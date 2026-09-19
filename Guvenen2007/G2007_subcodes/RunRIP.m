% RunRIP.m
% Stage (b) of the Guvenen (2007 AER) replication: solve the RIP model (no learning)
% through the FHorz Case1 machinery, as a plumbing rehearsal for the HIP baseline.
% Requires GPU. Run from baseline/.
%
% Exogenous states are a 3-dim z-block with age-dependent z_grid_J/pi_z_J:
%   z1 = alpha (fixed effect): identity transition at every age
%   z2 = z (AR(1) from z_0=0): discretizeLifeCycleAR1_KFTT over the 40 working ages
%   z3 = e (transitory): iid while working
% From the last working period on, ALL transitions become identity ("freeze"), so the
% pension Phi(Ytilde_T) can be computed inside the ReturnFn from the frozen states.
% This kron-assembled, age-dependent, partly-degenerate pi_z_J is exactly the structure
% the HIP belief chain needs in stage (c).
%
% Checks are inline (pre-solve chain checks, post-solve economics checks), printed
% PASS/FAIL to the diary. Do not silently retune anything: a FAIL is a report.

clearvars -except doPart
if exist('RunRIP_diary.txt','file'); delete('RunRIP_diary.txt'); end % diary appends; start fresh each run
if ~exist('SavedOutput/Graphs','dir'); mkdir('SavedOutput/Graphs'); end
diary('RunRIP_diary.txt');
fprintf('=== RunRIP.m run %s ===\n',char(datetime('now')));
npass=0; nfail=0;

%% Parameters
% Demographics: enter 25, retire 65, die 95; period = 1 year
Jwork=40;                    % working periods (ages 25-64)
Jret=31;                     % retirement periods (ages 65-95)
N_j=Jwork+Jret;              % 71
agevec=25:95;

% Income process: Table 1 row (1), RIP, all individuals
Params.rho=0.988;
Params.sigma2_alpha=0.058;
Params.sigma2_eta=0.015;
Params.sigma2_eps=0.061;
Params.betabar=0.009;        % mean income growth (p.701)
Params.alphabar=1.5;         % normalization (p.701)

% Preferences etc: Table 3
Params.crra=2;
Params.Pb=0.96;              % bond price; r = 1/Pb - 1 = 4.16%
Params.r=1/Params.Pb-1;
Params.delta=0.964;          % RIP time-discount factor (p.700; targets wealth/income=4)
Params.Phibar=1/1.40;        % pension scaling (footnote 18)

% Age-dependent basics
Params.exper_j=[0:Jwork-1, (Jwork-1)*ones(1,Jret)]; % experience, 0 at age 25; frozen in retirement
Params.exper_T=Jwork-1;      % experience at the last working period (used by pension)
Params.workret_j=[ones(1,Jwork), zeros(1,Jret)];
Params.mewj=ones(1,N_j)/N_j; % uniform age weights (no mortality risk in the paper)
AgeWeightParamNames={'mewj'};
DiscountFactorParamNames={'delta'};

% HIP income-process parameters (Table 1 row (2)): used ONLY for the natural borrowing
% limit, which the paper computes under HIP public-prior information and then imposes
% identically in the RIP model "for comparability" (p.702).
hip.rho=0.821; hip.s2a=0.022; hip.s2b=0.00038; hip.s2eta=0.029; hip.s2eps=0.047;

%% Grid sizes
n_d=0; d_grid=[];
% n_a is set below once the asset grid is built (121 negative + 381 positive points, 0 shared)
n_alpha=7; n_zshock=15; n_eps=5;
n_z=[n_alpha,n_zshock,n_eps];
N_z=prod(n_z);

%% Exogenous state grids
% Discretizers run on CPU (parallel=0) so the kron assembly below stays on CPU;
% the toolkit converts grids/transitions to GPU internally.
discopts.parallel=0;
discopts.nSigmas=2.5; % paper truncates the income distribution at 2.5 sd (p.702)
% alpha: iid normal quadrature, mean alphabar, var sigma2_alpha (identity transition)
[alpha_grid,w_alpha]=discretizeIIDNormal_TanakaToda(Params.alphabar,sqrt(Params.sigma2_alpha),n_alpha,discopts);
% e: transitory, mean 0, var sigma2_eps
[eps_grid,w_eps]=discretizeIIDNormal_TanakaToda(0,sqrt(Params.sigma2_eps),n_eps,discopts);
% z: life-cycle AR(1), z_0=0, over the working ages
[zshock_grid_J,pi_zshock_J,jequaloneDistz]=discretizeLifeCycleAR1_KFTT(...
    zeros(1,Jwork),Params.rho*ones(1,Jwork),sqrt(Params.sigma2_eta)*ones(1,Jwork),n_zshock,Jwork,discopts);
alpha_grid=gather(alpha_grid); w_alpha=gather(w_alpha);
eps_grid=gather(eps_grid); w_eps=gather(w_eps);
zshock_grid_J=gather(zshock_grid_J); pi_zshock_J=gather(pi_zshock_J); jequaloneDistz=gather(jequaloneDistz);

% Stacked age-dependent grid [sum(n_z), N_j]: column j = [alpha_grid; zshock_grid(:,j); eps_grid]
% Retirement ages reuse the last working-age z grid (frozen states)
z_grid_J=zeros(sum(n_z),N_j);
for jj=1:N_j
    jz=min(jj,Jwork);
    z_grid_J(:,jj)=[alpha_grid; zshock_grid_J(:,jz); eps_grid];
end

% Age-dependent transition [N_z,N_z,N_j-1]; slice jj = transition from period jj to jj+1.
% Convention: first z variable changes fastest, so joint pi = kron(pi_e, kron(pi_z, pi_alpha)).
pi_eps_iid=ones(n_eps,1)*w_eps'; % iid: every row is the e distribution
pi_z_J=zeros(N_z,N_z,N_j-1);
for jj=1:N_j-1
    if jj<=Jwork-1 % working-to-working (incl. into the last working period)
        pi_z_J(:,:,jj)=kron(pi_eps_iid,kron(pi_zshock_J(:,:,jj),eye(n_alpha)));
    else % last-working-to-retirement and within retirement: freeze everything
        pi_z_J(:,:,jj)=eye(N_z);
    end
end

% Joint grid at the last working age (for Ybar_T and checks); z1 fastest
jointgrid_T=[repmat(alpha_grid,n_zshock*n_eps,1), ...
    repmat(kron(zshock_grid_J(:,Jwork),ones(n_alpha,1)),n_eps,1), ...
    kron(eps_grid,ones(n_alpha*n_zshock,1))];

%% Ybar_T: mean income at the last working period (pension denominator)
% Push the initial z-distribution forward through the chain to age Jwork
pjoint=kron(w_eps,kron(jequaloneDistz,w_alpha))'; % 1-by-N_z, ordering matches kron above
for jj=1:Jwork-1
    pjoint=pjoint*pi_z_J(:,:,jj);
end
ygrid_T=exp(jointgrid_T(:,1)+Params.betabar*Params.exper_T+jointgrid_T(:,2)+jointgrid_T(:,3));
Params.Ybar_T=pjoint*ygrid_T;
fprintf('Ybar_T (mean income at last working age) = %.4f\n',Params.Ybar_T);

%% Natural borrowing limit (common, age-specific; HIP public-prior information, p.702)
% Worst-case income path: alpha, beta at -2.5sd (public prior: full sigma2_beta),
% every eta and eps innovation at -2.5sd (paper truncates income shocks at 2.5sd).
ymin=zeros(1,N_j);
zmin=0;
for jj=1:Jwork
    zmin=hip.rho*zmin-2.5*sqrt(hip.s2eta);
    ymin(jj)=exp((Params.alphabar-2.5*sqrt(hip.s2a))+(Params.betabar-2.5*sqrt(hip.s2b))*Params.exper_j(jj)+zmin-2.5*sqrt(hip.s2eps));
end
% Retirement: worst-case pension (bottom bracket of the replacement rule)
ymin(Jwork+1:N_j)=Params.Phibar*0.9*ymin(Jwork); % Ytilde_min=ymin(T)/Ybar_T, times Ybar_T cancels
% M(j) = PV at start of j of the minimum income stream from j onward
M=zeros(1,N_j+1);
for jj=N_j:-1:1
    M(jj)=ymin(jj)+Params.Pb*M(jj+1);
end
% Limit on the aprime chosen at age j: aprime >= -M(j+1); last period aprime>=0
Params.Wbar_j=[-M(2:N_j),0];
fprintf('Natural borrowing limit: min(Wbar)=%.3f (age %d), Wbar at last working age=%.3f\n',...
    min(Params.Wbar_j),agevec(find(Params.Wbar_j==min(Params.Wbar_j),1)),Params.Wbar_j(Jwork));

%% Asset grid: negative segment (down to the loosest borrowing limit) + curved positive segment
amin=min(Params.Wbar_j);
amax=60*Params.Ybar_T;
a_grid=unique([linspace(amin,0,97)'; (amax*(linspace(0,1,279).^3))']);
n_a=length(a_grid); % 375 (0 appears in both pieces)
fprintf('Asset grid: %d points on [%.2f, %.2f]\n',n_a,amin,amax);

%% Pre-solve checks on the exogenous-state objects
fprintf('\n--- Pre-solve checks ---\n');
% (i) every transition row sums to 1
maxrowdev=0;
for jj=1:N_j-1
    maxrowdev=max(maxrowdev,max(abs(sum(pi_z_J(:,:,jj),2)-1)));
end
fprintf('max |row sum - 1| over all pi_z_J slices = %.3e (tol 1e-12)\n',maxrowdev);
if maxrowdev<1e-12; npass=npass+1; else; nfail=nfail+1; fprintf('FAIL: pi_z_J row sums\n'); end
% (ii) retirement slices are exactly identity
maxidev=0;
for jj=Jwork:N_j-1
    maxidev=max(maxidev,max(abs(pi_z_J(:,:,jj)-eye(N_z)),[],'all'));
end
fprintf('max |retirement slice - identity| = %.3e (tol 0)\n',maxidev);
if maxidev==0; npass=npass+1; else; nfail=nfail+1; fprintf('FAIL: retirement freeze\n'); end
% (iii) chain-implied var(z_j) vs analytic sigma2_eta*(1-rho^2j)/(1-rho^2), and
% (iv) chain-implied var(log y_j) vs analytic (also validates the kron ordering)
pj=kron(w_eps,kron(jequaloneDistz,w_alpha))';
maxrelerr_vz=0; maxrelerr_vy=0;
for jj=1:Jwork
    if jj>1; pj=pj*pi_z_J(:,:,jj-1); end
    jointgrid_j=[repmat(alpha_grid,n_zshock*n_eps,1), ...
        repmat(kron(zshock_grid_J(:,jj),ones(n_alpha,1)),n_eps,1), ...
        kron(eps_grid,ones(n_alpha*n_zshock,1))];
    vz_chain=pj*(jointgrid_j(:,2).^2)-(pj*jointgrid_j(:,2))^2;
    vz_analytic=Params.sigma2_eta*(1-Params.rho^(2*jj))/(1-Params.rho^2);
    maxrelerr_vz=max(maxrelerr_vz,abs(vz_chain/vz_analytic-1));
    logy_j=jointgrid_j(:,1)+Params.betabar*Params.exper_j(jj)+jointgrid_j(:,2)+jointgrid_j(:,3);
    vy_chain=pj*(logy_j.^2)-(pj*logy_j)^2;
    vy_analytic=Params.sigma2_alpha+Params.sigma2_eps+vz_analytic;
    maxrelerr_vy=max(maxrelerr_vy,abs(vy_chain/vy_analytic-1));
end
fprintf('max rel err chain var(z_j) vs analytic over working ages = %.3e (tol 2e-2)\n',maxrelerr_vz);
fprintf('max rel err chain var(log y_j) vs analytic               = %.3e (tol 2e-2)\n',maxrelerr_vy);
if maxrelerr_vz<0.02; npass=npass+1; else; nfail=nfail+1; fprintf('FAIL: var(z_j) discretization\n'); end
if maxrelerr_vy<0.02; npass=npass+1; else; nfail=nfail+1; fprintf('FAIL: var(log y_j) (check kron ordering!)\n'); end
% (v) ReturnFn spot checks vs hand-computed values (CPU, scalar calls)
c1=0+exp(1.5+0+0)-Params.Pb*1; % working, alpha=1.5, exper=0, z=e=0, a=0, aprime=1
F1=ReturnFn_RIP(1,0,1.5,0,0,1,0,-10,Params.betabar,Params.Pb,Params.crra,Params.Phibar,Params.Ybar_T,Params.exper_T);
err1=abs(F1-(c1^(1-Params.crra))/(1-Params.crra));
Yt2=exp(1.5+Params.betabar*Params.exper_T)/Params.Ybar_T; % retired, z=e=0
if Yt2<0.3; Ph2=0.9*Yt2; elseif Yt2<=2; Ph2=0.27+0.32*(Yt2-0.3); elseif Yt2<=4.1; Ph2=0.81+0.15*(Yt2-2); else; Ph2=1.1; end
c2=2+Params.Phibar*Ph2*Params.Ybar_T-Params.Pb*1;
F2=ReturnFn_RIP(1,2,1.5,0,0,0,Params.exper_T,-10,Params.betabar,Params.Pb,Params.crra,Params.Phibar,Params.Ybar_T,Params.exper_T);
err2=abs(F2-(c2^(1-Params.crra))/(1-Params.crra));
F3=ReturnFn_RIP(-11,0,1.5,0,0,1,0,-10,Params.betabar,Params.Pb,Params.crra,Params.Phibar,Params.Ybar_T,Params.exper_T); % violates limit
fprintf('ReturnFn spot checks: working err=%.3e, retired err=%.3e, limit-violation F=%.1f (want -Inf)\n',err1,err2,F3);
if err1<1e-12 && err2<1e-12 && F3==-Inf; npass=npass+1; else; nfail=nfail+1; fprintf('FAIL: ReturnFn spot checks\n'); end
% (vi) Wbar series sensible: nondecreasing toward 0 late in life, 0 at the end
if all(diff(Params.Wbar_j(Jwork:end))>=-1e-12) && Params.Wbar_j(N_j)==0
    npass=npass+1; fprintf('Wbar series: retirement portion nondecreasing to 0: PASS\n');
else
    nfail=nfail+1; fprintf('FAIL: Wbar series shape\n');
end

%% Setup for solving
ReturnFn=@(aprime,a,alpha,z,e,workret_j,exper_j,Wbar_j,betabar,Pb,crra,Phibar,Ybar_T,exper_T) ...
    ReturnFn_RIP(aprime,a,alpha,z,e,workret_j,exper_j,Wbar_j,betabar,Pb,crra,Phibar,Ybar_T,exper_T);
vfoptions.verbose=0; % solve takes ~6s; per-age progress lines just clutter the diary
vfoptions.divideandconquer=1; % DC+GI: interpolate aprime between nodes
vfoptions.gridinterplayer=1;  % interpolate aprime between asset nodes (smooth policy)
vfoptions.ngridinterp=10;
simoptions.gridinterplayer=1; % MUST match vfoptions for StationaryDist/LifeCycleProfiles
simoptions.ngridinterp=10;

% Agent distribution at j=1: zero assets; alpha ~ its weights; z_1 ~ jequaloneDistz; e ~ its weights
[~,a0ind]=min(abs(a_grid)); % index of a=0
jequaloneDist=zeros([n_a,n_z],'gpuArray');
jequaloneDist(a0ind,:,:,:)=reshape(kron(w_eps,kron(jequaloneDistz,w_alpha)),[1,n_z]);

%% FnsToEvaluate (needed by the delta calibration loop as well as the final profiles)
FnsToEvaluate.cons=@(aprime,a,alpha,z,e,workret_j,exper_j,betabar,Pb,Phibar,Ybar_T,exper_T) ...
    a-Pb*aprime+workret_j*exp(alpha+betabar*exper_j+z+e)+(1-workret_j)*Phibar*Ybar_T*( ...
    (exp(alpha+betabar*exper_T+z+e)/Ybar_T<0.3)*(0.9*exp(alpha+betabar*exper_T+z+e)/Ybar_T) ...
    +(exp(alpha+betabar*exper_T+z+e)/Ybar_T>=0.3)*(exp(alpha+betabar*exper_T+z+e)/Ybar_T<=2)*(0.27+0.32*(exp(alpha+betabar*exper_T+z+e)/Ybar_T-0.3)) ...
    +(exp(alpha+betabar*exper_T+z+e)/Ybar_T>2)*(exp(alpha+betabar*exper_T+z+e)/Ybar_T<=4.1)*(0.81+0.15*(exp(alpha+betabar*exper_T+z+e)/Ybar_T-2)) ...
    +(exp(alpha+betabar*exper_T+z+e)/Ybar_T>4.1)*1.1);
FnsToEvaluate.income=@(aprime,a,alpha,z,e,workret_j,exper_j,betabar,Phibar,Ybar_T,exper_T) ...
    workret_j*exp(alpha+betabar*exper_j+z+e)+(1-workret_j)*Phibar*Ybar_T*( ...
    (exp(alpha+betabar*exper_T+z+e)/Ybar_T<0.3)*(0.9*exp(alpha+betabar*exper_T+z+e)/Ybar_T) ...
    +(exp(alpha+betabar*exper_T+z+e)/Ybar_T>=0.3)*(exp(alpha+betabar*exper_T+z+e)/Ybar_T<=2)*(0.27+0.32*(exp(alpha+betabar*exper_T+z+e)/Ybar_T-0.3)) ...
    +(exp(alpha+betabar*exper_T+z+e)/Ybar_T>2)*(exp(alpha+betabar*exper_T+z+e)/Ybar_T<=4.1)*(0.81+0.15*(exp(alpha+betabar*exper_T+z+e)/Ybar_T-2)) ...
    +(exp(alpha+betabar*exper_T+z+e)/Ybar_T>4.1)*1.1);
FnsToEvaluate.assets=@(aprime,a,alpha,z,e) a;
% log(max(c,1e-12)): FnsToEvaluate are evaluated on the WHOLE grid, including
% infeasible states with zero distribution mass where c<=0; GPU arrayfun errors on
% log of a negative. The guard never binds at any state with positive mass.
FnsToEvaluate.logc=@(aprime,a,alpha,z,e,workret_j,exper_j,betabar,Pb,Phibar,Ybar_T,exper_T) ...
    log(max(1e-12,a-Pb*aprime+workret_j*exp(alpha+betabar*exper_j+z+e)+(1-workret_j)*Phibar*Ybar_T*( ...
    (exp(alpha+betabar*exper_T+z+e)/Ybar_T<0.3)*(0.9*exp(alpha+betabar*exper_T+z+e)/Ybar_T) ...
    +(exp(alpha+betabar*exper_T+z+e)/Ybar_T>=0.3)*(exp(alpha+betabar*exper_T+z+e)/Ybar_T<=2)*(0.27+0.32*(exp(alpha+betabar*exper_T+z+e)/Ybar_T-0.3)) ...
    +(exp(alpha+betabar*exper_T+z+e)/Ybar_T>2)*(exp(alpha+betabar*exper_T+z+e)/Ybar_T<=4.1)*(0.81+0.15*(exp(alpha+betabar*exper_T+z+e)/Ybar_T-2)) ...
    +(exp(alpha+betabar*exper_T+z+e)/Ybar_T>4.1)*1.1)));
FnsToEvaluate.logy=@(aprime,a,alpha,z,e,workret_j,exper_j,betabar,exper_T) ...
    workret_j*(alpha+betabar*exper_j+z+e)+(1-workret_j)*(alpha+betabar*exper_T+z+e); % log labor income (frozen in retirement; use working ages only)


%% Calibrate delta to the wealth/income target of 4 (the paper's own procedure, p.700;
% the paper's resulting value is 0.964 -- we report ours and record any gap in NOTES.md).
% W/Y here uses total income (labor + pension + capital income r*a), matching the
% aggregate SCF/NIPA-style ratios of Budria Rodriguez et al. (2002); the labor+pension
% definition is also reported each iteration.
do_deltacalib=1; % =0 to just use Params.delta as set above
if do_deltacalib==1
    fprintf('\n--- Calibrating delta to W/Y=4 ---\n');
    FnsCalib.assets=FnsToEvaluate.assets;
    FnsCalib.income=FnsToEvaluate.income;
    dlo=0.90; dhi=0.978;
    for calibit=1:12
        Params.delta=(dlo+dhi)/2;
        [~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
        StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z_J,Params,simoptions);
        ACScalib=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsCalib,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,simoptions);
        meana_c=gather(ACScalib.assets.Mean); meany_c=gather(ACScalib.income.Mean);
        WYtot=sum(Params.mewj.*meana_c)/sum(Params.mewj.*(meany_c+Params.r*meana_c));
        WYlab=sum(Params.mewj.*meana_c)/sum(Params.mewj.*meany_c);
        fprintf('  delta=%.5f: W/Y total-income=%.3f (labor+pension def: %.3f)\n',Params.delta,WYtot,WYlab);
        if WYtot>4; dhi=Params.delta; else; dlo=Params.delta; end
    end
    Params.delta=(dlo+dhi)/2;
    fprintf('calibrated delta = %.4f (paper reports 0.964)\n',Params.delta);
end

%% Solve the value function (final)
fprintf('\n--- ValueFnIter ---\n');
tic;
[V,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
vftime=toc;
fprintf('ValueFnIter time: %.1f seconds\n',vftime);

% V monotone (nondecreasing) in assets at a few (age,z) slices.
% Deep-borrowing states at older ages can be infeasible (V=-Inf, zero dist mass);
% compare only the finite segment (diff of two -Inf is NaN, which would fail >=).
Vtest=reshape(gather(V),[n_a,N_z,N_j]);
mono=true;
for jj=[1,20,Jwork,Jwork+1,N_j]
    for zz=[1,ceil(N_z/2),N_z]
        vslice=Vtest(:,zz,jj);
        vslice=vslice(isfinite(vslice));
        mono=mono && all(diff(vslice)>=-1e-9);
    end
end
if mono; npass=npass+1; fprintf('V nondecreasing in assets: PASS\n');
else; nfail=nfail+1; fprintf('FAIL: V not monotone in assets\n'); end

%% Agent distribution (jequaloneDist built above)
fprintf('\n--- StationaryDist ---\n');
tic;
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z_J,Params,simoptions);
disttime=toc;
fprintf('StationaryDist time: %.1f seconds\n',disttime);

%% Life-cycle profiles
fprintf('\n--- LifeCycleProfiles ---\n');
tic;
AgeConditionalStats=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,simoptions);
proftime=toc;
fprintf('LifeCycleProfiles time: %.1f seconds\n',proftime);

%% Post-solve economics checks
fprintf('\n--- Post-solve checks ---\n');
meanc=AgeConditionalStats.cons.Mean;
meany=AgeConditionalStats.income.Mean;
meana=AgeConditionalStats.assets.Mean;
varlogc=AgeConditionalStats.logc.Variance;
varlogy=AgeConditionalStats.logy.Variance;

% (a) profile shapes. NB: this model has NO whole-life consumption hump: with
% delta*(1+r)>=1, no mortality risk and no bequest, mean consumption drifts up at the
% Euler rate (delta*(1+r))^(1/crra)-1 through retirement (first run: 0.2152%%/yr vs
% 0.2081%%/yr predicted). The paper only plots ages 25-65. So check instead:
% (a1) working-life consumption growth c(65)/c(25) near Fig 11's RIP value ~1.43
cgrowth=meanc(Jwork+1)/meanc(1);
fprintf('mean consumption: c(25)=%.3f, c(65)=%.3f, c(95)=%.3f; c(65)/c(25)=%.3f (paper Fig 11 RIP: ~1.43)\n',...
    meanc(1),meanc(Jwork+1),meanc(end),cgrowth);
if cgrowth>1.2 && cgrowth<1.75; npass=npass+1; else; nfail=nfail+1; fprintf('FAIL: working-life consumption growth off Fig 11 range\n'); end
% (a2) retirement consumption drift matches the unconstrained Euler equation
gret=mean(meanc(Jwork+5:N_j-1)./meanc(Jwork+4:N_j-2))-1;
geuler=(Params.delta*(1+Params.r))^(1/Params.crra)-1;
fprintf('retirement mean-c growth %.4f%%/yr vs Euler prediction %.4f%%/yr\n',100*gret,100*geuler);
if abs(gret-geuler)<0.001; npass=npass+1; else; nfail=nfail+1; fprintf('FAIL: retirement Euler drift\n'); end
% (a3) mean assets peak at retirement (life-cycle triangle)
[apk,iapk]=max(meana);
fprintf('mean assets: a(25)=%.2f, peak %.1f at age %d, a(95)=%.2f\n',meana(1),apk,agevec(iapk),meana(end));
if agevec(iapk)>=60 && agevec(iapk)<=70; npass=npass+1; else; nfail=nfail+1; fprintf('FAIL: asset peak not at retirement\n'); end
% (b) wealth/income ratio at the calibrated delta
WYtot=sum(Params.mewj.*meana)/sum(Params.mewj.*(meany+Params.r*meana));
WYlab=sum(Params.mewj.*meana)/sum(Params.mewj.*meany);
fprintf('aggregate W/Y: total-income def = %.2f (target 4), labor+pension def = %.2f; delta=%.4f (paper: 0.964)\n',WYtot,WYlab,Params.delta);
if WYtot>3.7 && WYtot<4.3; npass=npass+1; else; nfail=nfail+1; fprintf('FAIL: W/Y not at target after calibration\n'); end
% (c) var(log y) at ages 30/45/64 vs analytic (distribution-based, so also checks the dist)
for aa=[30 45 64]
    jj=aa-24;
    van=Params.sigma2_alpha+Params.sigma2_eps+Params.sigma2_eta*(1-Params.rho^(2*jj))/(1-Params.rho^2);
    fprintf('var(log y) age %d: model %.4f, analytic %.4f\n',aa,varlogy(jj),van);
end
vlyerr=abs(varlogy(45-24)/(Params.sigma2_alpha+Params.sigma2_eps+Params.sigma2_eta*(1-Params.rho^(2*21))/(1-Params.rho^2))-1);
if vlyerr<0.03; npass=npass+1; else; nfail=nfail+1; fprintf('FAIL: var(log y) from dist vs analytic\n'); end
% (d) rise in var(log c) over working life (Fig 7: RIP model rises ~0.26, 25 to 65)
dvarlogc=varlogc(Jwork+1)-varlogc(1);
fprintf('var(log c): %.3f at 25, %.3f at 65, rise = %.3f (paper Fig 7: ~0.26 rise)\n',varlogc(1),varlogc(Jwork+1),dvarlogc);
if dvarlogc>0.10 && dvarlogc<0.45; npass=npass+1; else; nfail=nfail+1; fprintf('FAIL: var(log c) rise off Fig 7 range\n'); end
% (e) no mass piling at the grid edges
topmass=gather(sum(StationaryDist(end-2:end,:,:,:,:),'all'));
botmass=gather(sum(StationaryDist(1:3,:,:,:,:),'all'));
fprintf('mass in top 3 asset points = %.3e (tol 1e-5); bottom 3 points = %.3e (info only)\n',topmass,botmass);
if topmass<1e-5; npass=npass+1; else; nfail=nfail+1; fprintf('FAIL: mass at top of asset grid, raise amax\n'); end

%% Figures
fig=figure(1); clf;
plot(agevec,meanc,'-','LineWidth',1.5); hold on; plot(agevec,meany,'--','LineWidth',1.5);
xlabel('age'); legend('mean consumption','mean income','Location','best'); title('RIP model: consumption and income');
saveas(fig,'SavedOutput/Graphs/RunRIP_fig1_meancy.png');
fig=figure(2); clf;
plot(agevec,meana,'-','LineWidth',1.5); xlabel('age'); title('RIP model: mean assets');
saveas(fig,'SavedOutput/Graphs/RunRIP_fig2_meanassets.png');
fig=figure(3); clf;
plot(agevec,varlogc,'-','LineWidth',1.5); hold on; plot(agevec(1:Jwork),varlogy(1:Jwork),'--','LineWidth',1.5);
xlabel('age'); legend('var(log c)','var(log y), working ages','Location','northwest');
title('RIP model: inequality profiles (cf. Fig 7)');
saveas(fig,'SavedOutput/Graphs/RunRIP_fig3_varlogs.png');

%% Save and summarize
save('RunRIP_results.mat','V','Policy','StationaryDist','AgeConditionalStats','Params','a_grid','z_grid_J','pi_z_J','n_a','n_z','N_j','-v7.3');
fprintf('\n=== RunRIP summary: %d passed, %d failed; runtimes vf=%.1fs dist=%.1fs prof=%.1fs ===\n',npass,nfail,vftime,disttime,proftime);
diary off;
