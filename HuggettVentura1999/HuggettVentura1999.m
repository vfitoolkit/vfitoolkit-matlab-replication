% Huggett & Ventura (1999) - On the Distributional Effects of Social Security Reform
% Review of Economic Dynamics 2, 498-531.
%
% Phase 1 + Phase 3: no idiosyncratic shocks. Solve
%   (A) U.S. system (concave PIA + common medical transfer b_common)
%   (B) Boskin proposal with floor b_lower = 0 (annuity proportional to a2)
% and compute compensating variations by ability level (Fig 3b).
%
% Both economies use the experienceassetz framework:
%   U.S.  : experience asset = ebar (average past earnings); aprime(l,ebar,z)
%   Boskin: experience asset = a2   (accumulated SS taxes ); aprime(l,a2,z)

% Which of the paper's two wage specifications to solve (Section 3.2).
%   0 = "no idiosyncratic shock": (gamma,sigma_eps^2)=(1.0,0.0), sigma_y1^2=0.376.
%       Ability is drawn once at birth and never moves -> PERMANENT TYPES.
%       This is the Table III / Fig 3b economy.
%   1 = "idiosyncratic shock":    (gamma,sigma_eps^2)=(0.985,0.015), sigma_y1^2=0.27.
%       z is a genuine AR(1) in each period of life -> an exogenous MARKOV STATE.
%       This is the Table IV / Fig 3a economy.
% The .mat outputs are tagged with the spec so the two specifications' results
% coexist -- which matters, because doPart(3) loads the .mat files and must not
% pick up the other specification's. The DIARY is a single name and overwrites.
idioShocks=1;
% specstr is just a filename tag on the .mat outputs, 'noidio' or 'idio'.

if idioShocks==0
    specstr='noidio';
else
    specstr='idio';
end

% Record the full run (incl. any error) to HV1999_Diary.txt, fresh each
% run (diary otherwise appends; if a run errors mid-way the diary stays on
% and captures the error, and the next run's header cleans up and restarts)
diary off
if exist('./HV1999_Diary.txt','file'), delete('./HV1999_Diary.txt'); end
diary ./HV1999_Diary.txt

% Subcodes (ReturnFns, aprimeFns, the US data build) live in ./HV1999subcodes/;
% MATLAB does not search subfolders of the working directory, so put it on the
% path. Note this does NOT change cwd, so the './SavedOutput/...' paths below
% still resolve to here.
addpath('./HV1999subcodes');

% To run on server, we need to tell Matlab where to find VFI Toolkit
addpath(genpath('./VFIToolkit-matlab/'));

% Output folders: .mat into ./SavedOutput/, figures into ./SavedOutput/Graphs/,
% auto-written latex table fragments into ./SavedOutput/LatexInputs/.
if ~exist('./SavedOutput','dir'); mkdir('./SavedOutput'); end
if ~exist('./SavedOutput/Graphs','dir'); mkdir('./SavedOutput/Graphs'); end
if ~exist('./SavedOutput/LatexInputs','dir'); mkdir('./SavedOutput/LatexInputs'); end

% Run one stage per cluster job; each stage saves its own .mat so a later
% stage can be run on its own.
doPart=[0,0,0,0,0,0,1];
% doPart(1): U.S. system steady state                      -> Table III/IV row 1
% doPart(2): Boskin proposal, all floor levels             -> Table III/IV rows 2-5
% doPart(3): compensating variations, Boskin               -> Fig 3a/3b
% doPart(4): Boskin w/medical transfer, all floor levels   -> Table VI
% doPart(5): compensating variations, Boskin w/medical     -> Fig 4a/4b
% doPart(6): age profiles of hours and consumption         -> Fig 5, Fig 6
% doPart(7): all figures, from the saved .mat files        -> Figs 1-6
%
% Part 1 [idio] is already on disk from the 2026-10-15 run (Table IV: K/Y 3.043
% vs 3.06, r 5.83 vs 5.8, income Gini 0.430 vs 0.43), so later parts run off
% HV1999_USsystem_idio.mat via the guards below.
%
% Parts 2 and 4 each solve one GE per floor level, so this is 8 Boskin solves per
% spec at roughly 15-20 min each. Part 6 needs no new solve (it reads the saved
% Policy and StationaryDist), and part 7 needs no model output at all beyond the
% .mat files, so both are cheap and safe to re-run on their own.

% Floor benefit levels, as multiples of output per person (HV1999 p.514: "a
% number of values for the floor benefit level varying from zero times output
% per person to 0.35 times output per person"). Tables III and IV report exactly
% these four. b_lower itself is a GE unknown, since Y is endogenous.
Params.bfloorfracs=[0,0.15,0.25,0.35];

% Which of those floor levels doPart(2) and doPart(4) actually SOLVE in this job.
% Params.bfloorfracs above still defines the table rows and figure legends, so
% stages 3, 5 and 7 expect all four .mat files to exist; this only selects which
% ones get computed now.
%
% Why this exists: on 2026-10-16 levels 1-3 each solved cleanly and level 4 then
% OOMed, on grids that a workspace diff confirmed were byte-identical across all
% four. Nothing about the problem grew. What differs is the allocator state
% inherited from three prior solves -- gpuDevice's AvailableMemory is CUDA's
% total free bytes, not the largest contiguous block, so a 4.88 GB request can
% fail with 40 GB nominally free. One level per job gives each solve a fresh
% process, which is the cleanest way to rule that out.
% One per stage, because the two stages are generally at different points: the
% plain Boskin levels may be partly solved while the medical variant has none.
doFloorsB=[1,1,1,1];   % doPart(2), plain Boskin
doFloorsM=[1,1,1,1];   % doPart(4), Boskin w/medical

% CheckTopOfGrid threshold: mass above this in the top gridpoints of an
% endogenous grid triggers a warning that the ceiling is binding.
checkTol=1e-5;

Params.headlessFigures=0; % set to 1 to render figures offscreen (safe for headless MATLAB); PNGs still saved
if Params.headlessFigures==1
    set(0,'DefaultFigureVisible','off');
end

%% Sizes
n_d=75;            % labor grid l in [0,1), spacing 0.0133.
                   % Cut from 101 on 2026-10-16 for the idioShocks==1 Boskin
                   % stage. The EV block of the ExpAssetz DC1_GI1 raw holds
                   % several [N_d2*N_a1,N_a2,N_z,N_z] copies, so this scales
                   % linearly in n_d: with n_a_B=[251,50] it gives 3.32 GB per
                   % array, between the 2.24 GB that solves cleanly (n_d=101,
                   % N_a2=25) and the 4.47 GB that OOMed (n_d=101, N_a2=50).
                   % Untested middle ground, run to find out.
                   % Accuracy is not the binding consideration: 101 was
                   % adopted on 2026-10-06 purely to test
                   % whether the hours gap was grid-driven. It is not -- every
                   % aggregate was identical to 4 significant figures at 51 vs
                   % 101 (K/Y 2.850 vs 2.851, hours 0.351 both, CEV at ability
                   % 21 +12.19 both) -- so 51 is available as a free fallback
                   % once the OOM is understood. NB the paper does not grid the
                   % labour choice at all: its appendix (p.528-529) uses
                   % golden-section search over l.
n_a=[251,25];      % [assets, ebar]. Was [101,4]: 4 ebar points cannot carry an
                   % accumulating average (see the ebar grid construction below).
N_j=81;            % ages 20..100
% N_i and n_z are NOT set here: they depend on which wage specification is being
% solved and are set in the "Wage process (Section 3.2)" block below.
%   idioShocks==0 : N_i=21 permanent ability types, n_z=0 (no exogenous state)
%   idioShocks==1 : N_i=1,  n_z=21 (ability is an AR(1) state)

%% Parameters - HV1999 Table I
Params.agejshifter=19;
Params.J=N_j;
Params.agej=(1:1:N_j)';
Params.Jr=46;                 % retirement at real-life age 65 = period 46

% Preferences
Params.beta=1.011;
Params.sigma=2.0;
Params.nu=0.33;

% Production
Params.A=0.895944;
Params.alpha=0.36;
Params.delta=0.06;
Params.g=0.021;

% Demographics
Params.n=0.012;

% Government consumption
Params.GoverY=0.195;

%% US male survival probs sj and Hansen-style age-earnings profile kappa_j
HuggettVentura1999_USdata;
Params.sj=HV1999_sj;
Params.kappa_j=HV1999_ybar;
clear HV1999_sj HV1999_ybar

%% Age weights mewj (Section 2.1 footnote 9)
AgeWeightsParamNames={'mewj'};
Params.mewj=ones(N_j,1);
for jj=2:N_j
    Params.mewj(jj)=Params.mewj(jj-1)*Params.sj(jj-1)/(1+Params.n);
end
Params.mewj=Params.mewj/sum(Params.mewj);
Params.mewj=Params.mewj';

%% Wage process (Section 3.2). Two specifications, selected by idioShocks.
% Both are the same process, eqn (13):
%     y_j - kappa_j = gamma*(y_{j-1} - ybar_{j-1}) + eps_j,   eps ~ N(0,sigma_eps^2)
%     y_1 ~ N(ybar_1, sigma_y1^2),   z := y_j - kappa_j,   e(z,j) = exp(z + kappa_j)
% so z is the deviation of log endowment from its age mean, and productivity is
% exp(z) times the age profile. The two specs differ only in (gamma, sigma_eps^2).
%
% Discretization follows footnote 17 in BOTH cases: 21 evenly spaced values for z
% between -4*sigma_y1 and +4*sigma_y1, with probabilities from integrating the
% normal (i.e. Tauchen).
if idioShocks==0
    % gamma=1, sigma_eps=0: z never moves after birth. That is a permanent type,
    % not a shock, so it is carried as N_i types with NO z state at all -- which
    % drops the ability dimension from each type's value function and removes the
    % n_z-by-n_z expectation an identity pi_z would otherwise compute for nothing.
    Params.sigma_y1_sq=0.376;
    sigma_y1=sqrt(Params.sigma_y1_sq);
    N_i=21;
    n_z=0; z_grid=[]; pi_z=[];
    logability_grid=linspace(-4*sigma_y1,4*sigma_y1,N_i)';
    ability_edges=[-Inf;(logability_grid(1:end-1)+logability_grid(2:end))/2;Inf];
    Params.abilityi_dist=diff(normcdf(ability_edges,0,sigma_y1));
    Params.abilityi_dist=Params.abilityi_dist/sum(Params.abilityi_dist);
    PTypeDistParamNames={'abilityi_dist'};
    % LOG ability per type. The ReturnFn/aprimeFn take log productivity and
    % exponentiate internally, so this parameter slots into exactly the position
    % that z occupies as a state under idioShocks==1 -- one set of subcode files
    % serves both specs. NB deliberately NOT renormalized to E[exp(z)]=1: fn17
    % describes only the discretization and sigma_y1^2 is itself calibrated to a
    % wage Gini, which pins the spread rather than the mean.
    Params.logabilityi=logability_grid;
    Names_i=cell(1,N_i);
    for ii=1:N_i
        Names_i{ii}=sprintf('ability%02d',ii);
    end
else
    % gamma=0.985, sigma_eps^2=0.015, sigma_y1^2=0.27 (p.510, and Table II's
    % middle row). z moves every period, so it is a genuine exogenous state and
    % there are no permanent types: N_i=1 is a formality that keeps every command
    % call identical to the noidio branch.
    Params.sigma_y1_sq=0.27;
    Params.gamma_z=0.985;
    Params.sigma_eps_sq=0.015;
    sigma_y1=sqrt(Params.sigma_y1_sq);
    n_z=21;
    % Grid bounds are +/-4*sigma_y1 per fn17 -- NOT +/-4 stationary sd. The
    % process is non-stationary over the life cycle (var grows from 0.27 at birth
    % toward sigma_eps^2/(1-gamma^2)=0.503), so this grid is about 3.1 sd wide at
    % age 45. That is the paper's stated choice and is kept.
    z_grid=linspace(-4*sigma_y1,4*sigma_y1,n_z)';
    % Tauchen transitions for z' = gamma*z + eps, integrating the normal between
    % midpoints of the grid (open-ended at both ends so each row sums to 1).
    pi_z=zeros(n_z,n_z);
    zmid=(z_grid(1:end-1)+z_grid(2:end))/2;
    for zz=1:n_z
        edges=[-Inf;zmid;Inf];
        pi_z(zz,:)=diff(normcdf(edges,Params.gamma_z*z_grid(zz),sqrt(Params.sigma_eps_sq)))';
        pi_z(zz,:)=pi_z(zz,:)/sum(pi_z(zz,:));
    end
    % Age-1 distribution of z: N(0,sigma_y1^2) binned on the same grid.
    z_edges=[-Inf;zmid;Inf];
    jequaloneDistz=diff(normcdf(z_edges,0,sigma_y1));
    jequaloneDistz=jequaloneDistz/sum(jequaloneDistz);
    % No permanent heterogeneity in this spec.
    N_i=1;
    Names_i={'all'};
    Params.abilityi_dist=1;
    PTypeDistParamNames={'abilityi_dist'};
end

%% Grids (asset, ebar, labor)
d_grid=linspace(0,1-1e-3,n_d)';

% Asset grid: a_k = bar_a*(k/n_a(1))^2.35 (HV1999 appendix curvature)
% amax=40, lowered from 50. At the first solve the asset distribution had mean
% 1.619 and median 0.843, with p99=10.56 and p99.9=19.65 -- so 121 of the 251
% gridpoints sat above p99 holding under 1% of mass, while only 44 lay below the
% median. Lowering the ceiling redistributes points into the region where the
% agents actually are: 49 below the median and 64 below the mean at amax=40,
% truncating 2.0e-05 of mass (against 1.1e-05 already in the top ten points at
% amax=50). The paper does not pin this down -- its appendix gives only "an
% upper bound imposed on the asset space" with no number -- so this is our
% choice, not a replication target. Candidates measured from the solved
% distribution: amax=30 gives 55 below the median truncating 1.0e-04, amax=25
% gives 60 truncating 3.1e-04.
Params.amax=40;
k_idx=(0:n_a(1)-1)';
asset_grid=Params.amax*(k_idx/(n_a(1)-1)).^2.35;

% ebar grid. We work with ebar in "ratio to mean working-age earnings" units,
% so the grid is literal and w-free, and the PIA bend points are {0.20, 1.24}
% with the AIME cap at 2.47.
%
% Those four values are where the benefit formula kinks, so they are kept as
% exact nodes, but they are NOT a sufficient grid on their own: ebar is an
% accumulating average,
%     ebar' = [ebar*(agej-1) + min(e(z,j)*l/ybar_mean, 2.47)] / agej,
% which the experience-asset machinery interpolates onto this grid every
% period. With only the four kink values the whole accumulation path would be
% crushed onto them. So we span [0, 2.47] with 25 points, piecewise-uniform
% between the kinks, and weight the interior segment where most ebar mass sits.
% Point counts per segment include both endpoints, and the two shared endpoints
% are removed by unique(), giving 5+13+9-2 = 25.
% ebar is AIME expressed in multiples of AVERAGE EARNINGS, which is what puts
% the PIA bend points at 0.20 and 1.24 and the cap at 2.47 (p.513). Average
% earnings is an EQUILIBRIUM object here, because labour supply is endogenous:
% it is w*E[ability*kappa_j*l] per worker, not w*mean(kappa_j). Treating it as the
% latter made the normalizer 2.411x too large at the 2026-10-05 solve, which put
% the bend points and cap at an effective 0.48, 2.99 and 5.95 times mean
% earnings and compressed the Fig 3b welfare profile by about half. So
% ybar_mean_unwage is now a GE unknown (it was always meant to be -- see the
% "ybar_mean endogenous" comment on the Boskin GE block) pinned by
%     ybar_mean_unwage = Earnings/(w*workshare)
% i.e. average earnings per WORKER, divided by w so the ReturnFn's
% w*ybar_mean_unwage is average earnings itself. Per worker rather than per
% person follows U.S. practice, where the bend points are indexed to the
% average wage index; the two differ by 1/workshare, about 15% here.
% HV2000 avoided needing this by renormalizing E[z]=1, which only works because
% it has no labour choice and so mean earnings equal w*ybar exactly.
Params.workshare=sum(Params.mewj(Params.agej<Params.Jr));
% Initial guess only; the GE solves for it. mean(kappa_j) scaled by a plausible
% average labour supply, which is the right order (the solved value was ~0.45).
Params.ybar_mean_unwage=0.35*sum(Params.mewj((Params.agej<Params.Jr)).*Params.kappa_j(Params.agej<Params.Jr)') ...
    /Params.workshare;
ebar_grid=unique([linspace(0,0.20,5),linspace(0.20,1.24,13),linspace(1.24,2.47,9)])';
if numel(ebar_grid)~=n_a(2)
    error('ebar_grid has %d points but n_a(2)=%d',numel(ebar_grid),n_a(2))
end

a_grid=[asset_grid;ebar_grid];

%% experience asset setup
% ebar's law of motion depends on productivity. Under idioShocks==1 productivity
% is the exogenous STATE z, which is what "experienceassetz" means; under
% idioShocks==0 it is a per-type parameter and there is no z at all, so the plain
% "experienceasset" raw is the right one. Both pass log productivity into
% ebarprimeFn's 3rd slot.
if idioShocks==0
    vfoptions.experienceasset=1;
    simoptions.experienceasset=1;
    vfoptions.aprimeFn=@(l,ebar,agej,Jr,kappa_j,ybar_mean_unwage,logabilityi) ...
        HuggettVentura1999_ebarprimeFn(l,ebar,logabilityi,agej,Jr,kappa_j,ybar_mean_unwage);
else
    vfoptions.experienceassetz=1;
    simoptions.experienceassetz=1;
    % experienceassetz needs the z grid on simoptions too: the agent-dist and
    % stats codes rebuild z_gridvals_J themselves (via ExogShockSetup_FHorz)
    % rather than inheriting it from the solver, and error out without it.
    simoptions.z_grid=z_grid;
    vfoptions.aprimeFn=@(l,ebar,z,agej,Jr,kappa_j,ybar_mean_unwage) ...
        HuggettVentura1999_ebarprimeFn(l,ebar,z,agej,Jr,kappa_j,ybar_mean_unwage);
    % Loop the value fn solve over current z. Needed because the ExpAssetz
    % DC1_GI1 raw's terminal-age block leaves two [N_d2,n2long,N_a1,N_a2,N_z]
    % arrays (a1primeindexes, ReturnMatrix_ii) resident for the whole age loop
    % -- 6.71 GB each at these grid sizes -- which is what made the first-age
    % EV line OOM on 2026-10-09. lowmemory=1 drops the z dim from both (0.32 GB
    % each), freeing 12.8 GB. It does NOT shrink the EV block itself, which is
    % the remaining peak and scales as N_z^2.
    vfoptions.lowmemory=1;
end
simoptions.aprimeFn=vfoptions.aprimeFn;
simoptions.a_grid=a_grid;
simoptions.d_grid=d_grid;

% VFI solver options
vfoptions.divideandconquer=1;
vfoptions.gridinterplayer=1;
vfoptions.ngridinterp=30;
simoptions.gridinterplayer=vfoptions.gridinterplayer;
simoptions.ngridinterp=vfoptions.ngridinterp;

%% Initial GE prices
KoverL=(Params.alpha*Params.A/(0.06+Params.delta))^(1/(1-Params.alpha));
Params.r=Params.alpha*Params.A*KoverL^(Params.alpha-1)-Params.delta;
Params.w=(1-Params.alpha)*Params.A*KoverL^Params.alpha;
Params.tau=0.20;
Params.theta=0.10;
Params.b_common=0.1242*Params.A*KoverL^Params.alpha;  % 0.1242*Y_per_capita; Y refined in GE

%% Discount factor (HV1999 eqn 4: beta*(1+g)^{nu*(1-sigma)}*s_{j+1})
Params.beta_g=Params.beta*(1+Params.g)^(Params.nu*(1-Params.sigma));
DiscountFactorParamNames={'beta_g','sj'};

%% Return function, newborn distribution, and FnsToEvaluate
% The 5th argument of the underlying ReturnFn is LOG productivity: a per-type
% parameter when ability is permanent, the z state when it is not.
if idioShocks==0
    ReturnFn=@(l,aprime,a,ebar,logabilityi,agej,Jr,sigma,nu,g,r,w,tau,theta,kappa_j,b_common,ybar_mean_unwage) ...
        HuggettVentura1999_USsystem_ReturnFn(l,aprime,a,ebar,logabilityi,agej,Jr,sigma,nu,g,r,w,tau,theta,kappa_j,b_common,ybar_mean_unwage);
    % Every type born at (a=0, ebar=0); the spread across ability is the ptype
    % distribution, not part of this object.
    jequaloneDist=zeros(n_a);
    jequaloneDist(1,1)=1;
    FnsToEvaluate.K=@(l,aprime,a,ebar) a;
    FnsToEvaluate.L=@(l,aprime,a,ebar,agej,Jr,kappa_j,logabilityi) (agej<Jr)*exp(logabilityi)*kappa_j*l;
    FnsToEvaluate.Earnings=@(l,aprime,a,ebar,agej,Jr,kappa_j,w,logabilityi) (agej<Jr)*exp(logabilityi)*kappa_j*l*w;
    FnsToEvaluate.AccBeq=@(l,aprime,a,ebar,sj,r,tau,g) (1-sj)*aprime*(1+g)*(1+r*(1-tau));
    FnsToEvaluate.SSBenefits=@(l,aprime,a,ebar,agej,Jr,g,b_common,ybar_mean_unwage,w) ...
        (agej>=Jr)*( b_common + w*ybar_mean_unwage*( ...
              0.90*min(ebar,0.20) ...
            + 0.32*max(0,min(ebar,1.24)-0.20) ...
            + 0.15*max(0,ebar-1.24) ) / ((1+g)^(agej-Jr)) );
    FnsToEvaluate.Ebar=@(l,aprime,a,ebar) ebar;
    FnsToEvaluate.HoursWorked=@(l,aprime,a,ebar,agej,Jr) (agej<Jr)*l;
else
    ReturnFn=@(l,aprime,a,ebar,z,agej,Jr,sigma,nu,g,r,w,tau,theta,kappa_j,b_common,ybar_mean_unwage) ...
        HuggettVentura1999_USsystem_ReturnFn(l,aprime,a,ebar,z,agej,Jr,sigma,nu,g,r,w,tau,theta,kappa_j,b_common,ybar_mean_unwage);
    % Born at (a=0, ebar=0) with z drawn from the age-1 distribution.
    jequaloneDist=zeros([n_a,n_z]);
    jequaloneDist(1,1,:)=reshape(jequaloneDistz,[1,1,n_z]);
    FnsToEvaluate.K=@(l,aprime,a,ebar,z) a;
    FnsToEvaluate.L=@(l,aprime,a,ebar,z,agej,Jr,kappa_j) (agej<Jr)*exp(z)*kappa_j*l;
    FnsToEvaluate.Earnings=@(l,aprime,a,ebar,z,agej,Jr,kappa_j,w) (agej<Jr)*exp(z)*kappa_j*l*w;
    FnsToEvaluate.AccBeq=@(l,aprime,a,ebar,z,sj,r,tau,g) (1-sj)*aprime*(1+g)*(1+r*(1-tau));
    FnsToEvaluate.SSBenefits=@(l,aprime,a,ebar,z,agej,Jr,g,b_common,ybar_mean_unwage,w) ...
        (agej>=Jr)*( b_common + w*ybar_mean_unwage*( ...
              0.90*min(ebar,0.20) ...
            + 0.32*max(0,min(ebar,1.24)-0.20) ...
            + 0.15*max(0,ebar-1.24) ) / ((1+g)^(agej-Jr)) );
    FnsToEvaluate.Ebar=@(l,aprime,a,ebar,z) ebar;
    FnsToEvaluate.HoursWorked=@(l,aprime,a,ebar,z,agej,Jr) (agej<Jr)*l;
end

%% General equilibrium setup (HV1999 Section 2.3 conditions 1-6)
GEPriceParamNames={'r','w','tau','theta','b_common','ybar_mean_unwage'};

heteroagentoptions.verbose=1;

heteroagentoptions.intermediateEqns.Y=@(A,alpha,K,L) A*K^alpha*L^(1-alpha);

GeneralEqmEqns.CapitalMarket=@(r,A,alpha,K,L,delta) r-(alpha*A*K^(alpha-1)*L^(1-alpha)-delta);
GeneralEqmEqns.LaborMarket =@(w,A,alpha,K,L)       w-(1-alpha)*A*K^alpha*L^(-alpha);
% Govt budget (eqm cond 5, US system: TR=0): G = tau*(rK+wL) + T - 0
%   => tau*(rK+wL) + T - GoverY*Y = 0
GeneralEqmEqns.TaxRate    =@(tau,r,K,w,L,AccBeq,GoverY,Y) tau*(r*K+w*L)+AccBeq-GoverY*Y; % Section 2.3, eqm cond 5: G = tau*(rK+wL)+T-TR with T the accidental-bequest revenue "which the government taxes fully" (p.508). No T-hat lump-sum back to agents (cf. agent budget eqn 5 on p.504, which has no T term), so bequests just offset tax-financed G.
% SS budget (eqm cond 6): theta*w*L = aggregate benefits
GeneralEqmEqns.SSBalance  =@(theta,w,L,SSBenefits) theta*w*L-SSBenefits;
% Common medical transfer indexed to Y (Section 3.3 with US system)
GeneralEqmEqns.bCommon    =@(b_common,Y) b_common-0.1242*Y;
% Average earnings per worker, the unit the PIA bend points are expressed in.
GeneralEqmEqns.YbarMean   =@(ybar_mean_unwage,Earnings,w,workshare) ybar_mean_unwage-Earnings/(w*workshare);

%% GE solver: Anderson acceleration of the shooting map (fminalgo=9).
% The default Nelder-Mead needs a few hundred evaluations on a problem this
% size; each evaluation here is a full VFI over 131,775 states x 12,801 choices
% per age (n_d=51 makes this ~6x HV2000's heaviest model), so the default is
% not affordable. Anderson solved HV2000's 6-price problem in ~10 iterations.
%
% howtoupdate rows are {GEeqnName, PriceName, add, factor}. With add=0 the
% shooting step is p <- p - factor*residual, which is a damped move toward the
% target ONLY IF the residual is written as "price minus target" (factor=1 is
% then a full replacement). All five conditions above have a residual that is
% increasing in their own price, so add=0 is right for each:
%   CapitalMarket, LaborMarket, bCommon are literally price-minus-target.
%   TaxRate residual is (r*K+w*L)*(tau - tau_implied) and SSBalance residual is
%   (w*L)*(theta - theta_implied), i.e. price-minus-target scaled by a tax base
%   that is O(1) in these units, so factor=0.5 is still a damped half-step.
% r is damped harder because the Aiyagari fixed point (higher r -> more saving
% -> more K -> lower MPK) oscillates undamped.
heteroagentoptions.fminalgo=9;
heteroagentoptions.fminalgo9.howtoupdate={ ...
    'CapitalMarket', 'r',        0, 0.3; ...
    'LaborMarket',   'w',        0, 0.5; ...
    'TaxRate',       'tau',      0, 0.5; ...
    'SSBalance',     'theta',    0, 0.5; ...
    'bCommon',       'b_common', 0, 0.5; ...
    'YbarMean', 'ybar_mean_unwage', 0, 0.5};
heteroagentoptions.anderson.memory=5;
heteroagentoptions.anderson.warmup=2;
heteroagentoptions.anderson.maxiter=200;

%% Solve U.S. system steady state
if doPart(1)==1
fprintf('\n=== HV1999 U.S. system [%s] ===\n',specstr);

% Signature is (...,N_j, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn,...): n_p
% comes before pi_z, and pi_z before the grids. Returns [p_eqm, GEcondns].
GEstagelabel='U.S. system';
% Memory report immediately before the GE call. The value-function iteration
% OOMed here on 2026-10-09 under idioShocks==1, inside
% ValueFnIter_FHorz_ExpAssetz_DC1_GI1_nod1_raw, where five arrays of shape
% [N_d2*N_a1, N_a2, N_z, N_z] are live at once (Vlower, Vupper, aprimeProbs,
% skipinterp, EV). whos shows what THIS workspace is holding; gpuDevice shows
% what is actually free on the card, which is the number that matters -- a
% workspace variable that is a gpuArray occupies the card, and HV2000 was once
% OOMed by exactly that (stale V/Policy/StationaryDist left from a prior stage).
fprintf('\n---- memory before the %s GE solve ----\n',GEstagelabel);
try
    gd=gpuDevice;
    fprintf('GPU %s: %.2f GB free of %.2f GB total\n', ...
        gd.Name,gd.AvailableMemory/1e9,gd.TotalMemory/1e9);
catch ME
    fprintf('gpuDevice() unavailable: %s\n',ME.message);
end
whos
fprintf('---- end memory report ----\n\n');
[p_eqm,GeneralEqmCondn]=HeteroAgentStationaryEqm_Case1_FHorz_PType( ...
    n_d,n_a,n_z,N_j,Names_i,[],pi_z,d_grid,a_grid,z_grid,jequaloneDist, ...
    ReturnFn,FnsToEvaluate,GeneralEqmEqns,Params,DiscountFactorParamNames, ...
    AgeWeightsParamNames,PTypeDistParamNames,GEPriceParamNames, ...
    heteroagentoptions,simoptions,vfoptions);

% Loop over GEPriceParamNames rather than hand-listing: a hand-listed block
% silently goes stale the moment a price is added. It did -- ybar_mean_unwage
% was added to the price list on 2026-10-05 without a write-back line, so the
% 2026-10-05 solve re-solved V, the distribution and every statistic at the
% INITIAL GUESS 0.3769 instead of the equilibrium 0.4518, a 17% error in the
% object that sets the PIA bend points. Every number from that run is void.
for pp=1:length(GEPriceParamNames)
    Params.(GEPriceParamNames{pp})=p_eqm.(GEPriceParamNames{pp});
end

% Final solve at eqm prices
[V_US,Policy_US]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
StationaryDist_US=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,Policy_US,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);
AggVars_US=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist_US,Policy_US,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);

%% Distributional statistics (HV1999 Table III/IV: Income Gini; Section 3.2: wage Gini)
% Kept OUT of the GE FnsToEvaluate so the GE loop does not pay for them.
% Income per HV1999 footnote 24 (p.516): "labor plus asset income before taxes,
% plus social security transfers" -- so earnings gross of tau and theta, plus
% r*a gross of tau, plus the benefit.
% Wage is the wage RATE w*ability*kappa_j, and sigma_y1^2=0.376 was calibrated so
% its Gini is 0.35 among under-65s (p.511), hence the conditional restriction.
FnsToEvaluate_stats=FnsToEvaluate;
if idioShocks==0
    FnsToEvaluate_stats.Income=@(l,aprime,a,ebar,agej,Jr,r,kappa_j,w,g,b_common,ybar_mean_unwage,logabilityi) ...
        (agej<Jr)*exp(logabilityi)*kappa_j*l*w + r*a ...
        + (agej>=Jr)*( b_common + w*ybar_mean_unwage*( ...
              0.90*min(ebar,0.20) + 0.32*max(0,min(ebar,1.24)-0.20) + 0.15*max(0,ebar-1.24) ) ...
              / ((1+g)^(agej-Jr)) );
    FnsToEvaluate_stats.Wage=@(l,aprime,a,ebar,agej,Jr,kappa_j,w,logabilityi) w*exp(logabilityi)*kappa_j;
else
    FnsToEvaluate_stats.Income=@(l,aprime,a,ebar,z,agej,Jr,r,kappa_j,w,g,b_common,ybar_mean_unwage) ...
        (agej<Jr)*exp(z)*kappa_j*l*w + r*a ...
        + (agej>=Jr)*( b_common + w*ybar_mean_unwage*( ...
              0.90*min(ebar,0.20) + 0.32*max(0,min(ebar,1.24)-0.20) + 0.15*max(0,ebar-1.24) ) ...
              / ((1+g)^(agej-Jr)) );
    FnsToEvaluate_stats.Wage=@(l,aprime,a,ebar,z,agej,Jr,kappa_j,w) w*exp(z)*kappa_j;
end
simoptions_stats=simoptions;
if idioShocks==0
    simoptions_stats.conditionalrestrictions.Under65=@(l,aprime,a,ebar,agej,Jr) (agej<Jr);
else
    simoptions_stats.conditionalrestrictions.Under65=@(l,aprime,a,ebar,z,agej,Jr) (agej<Jr);
end
AllStats_US=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_US,Policy_US,FnsToEvaluate_stats,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_stats);
% The income Gini target is spec-specific: Table III is the no-idiosyncratic
% economy, Table IV the idiosyncratic one. Reporting Table III's 0.40 against an
% idioShocks==1 run reads as a 7% overshoot when the run is in fact exact.
if idioShocks==0
    giniRef_US=0.40; tableRef='Table III';
else
    giniRef_US=0.43; tableRef='Table IV';
end
fprintf('\nDistributional checks (US system):\n');
fprintf('  Income Gini (all agents)   = %.3f   [HV1999 %s: %.2f]\n',gather(AllStats_US.Income.Gini),tableRef,giniRef_US);
% 0.35 is a construction target, not an outcome: the paper picks gamma_z so the
% under-65 wage Gini equals 0.35 at each sigma_eps_sq (Table II's gamma column).
% Our idio implementation lands at 0.370; accepted as a known residual.
fprintf('  Wage Gini (under 65)       = %.3f   [calibration target 0.35]\n',gather(AllStats_US.Under65.Wage.Gini));

%% CheckTopOfGrid: mass sitting at the top of each endogenous grid.
% StationaryDist_US is [n_a(1) assets, n_a(2) ebar, n_z, N_j], so dim 1 is the
% asset grid and dim 2 is ebar.
%
% The two dimensions mean different things here. Mass piled at the top of the
% ASSET grid is a numerical artifact: households want to accumulate past
% Params.amax and cannot, so the ceiling is distorting the solution and amax
% should be widened. Mass at the top of the EBAR grid is not -- its four points
% are the PIA bend points and the AIME cap {0, 0.20, 1.24, 2.47}, so the top
% point IS the cap, and mass there is the share of agents whose benefit is
% capped. That is an economic feature of the U.S. system, so it is reported for
% information and never warned on.
% Under PType, StationaryDist_US is a struct: one [n_a(1) x n_a(2) x N_j] array
% per permanent type, each summing to 1, plus the mixing weights in .ptweights.
% Pool over types (and over age) before taking the marginals.
pooled_US=zeros(n_a);
for ii=1:N_i
    di=gather(StationaryDist_US.(Names_i{ii}));
    pooled_US=pooled_US+StationaryDist_US.ptweights(ii)*sum(di,3:ndims(di)); % sum over age AND, under idioShocks==1, over z -- leaving the two endogenous dims
end
mass_by_a_US   =sum(pooled_US,2);
mass_by_ebar_US=sum(pooled_US,1);
CheckTopOfGrid_US.MassAtTop10gridptsAssets=gather(sum(mass_by_a_US(end-9:end)));
CheckTopOfGrid_US.MassAtTop5gridptsAssets =gather(sum(mass_by_a_US(end-4:end)));
CheckTopOfGrid_US.MassAtTop2gridptsAssets =gather(sum(mass_by_a_US(end-1:end)));
CheckTopOfGrid_US.MassAtAIMEcap           =gather(mass_by_ebar_US(end));
fprintf('\nCheckTopOfGrid (US system), asset ceiling amax=%g, threshold %.0e:\n',Params.amax,checkTol);
fprintf('  mass in top 10 asset gridpts : %.4e\n',CheckTopOfGrid_US.MassAtTop10gridptsAssets);
fprintf('  mass in top  5 asset gridpts : %.4e\n',CheckTopOfGrid_US.MassAtTop5gridptsAssets);
fprintf('  mass in top  2 asset gridpts : %.4e\n',CheckTopOfGrid_US.MassAtTop2gridptsAssets);
fprintf('  share at the AIME cap        : %.4f   (economic, not a grid limit)\n',CheckTopOfGrid_US.MassAtAIMEcap);
if CheckTopOfGrid_US.MassAtTop10gridptsAssets>checkTol
    warning('US system: %.4e mass in top 10 asset gridpts (>%.0e); widen Params.amax',CheckTopOfGrid_US.MassAtTop10gridptsAssets,checkTol);
end

%% Save Phase-1 output
save(sprintf('./SavedOutput/HV1999_USsystem_%s.mat',specstr),'Params','V_US','Policy_US','StationaryDist_US','AggVars_US','AllStats_US','p_eqm','GeneralEqmCondn','CheckTopOfGrid_US');

% Headline aggregates (Table III row 1)
Y_US=Params.A*AggVars_US.K.Mean^Params.alpha*AggVars_US.L.Mean^(1-Params.alpha);
KoverY_US=AggVars_US.K.Mean/Y_US;
fprintf('\nUS system [%s] headline aggregates:\n',specstr);
fprintf('  K/Y          = %.3f\n', KoverY_US);
fprintf('  L            = %.3f\n', AggVars_US.L.Mean);
fprintf('  HoursWorked  = %.3f\n', AggVars_US.HoursWorked.Mean / sum(Params.mewj(Params.agej<Params.Jr)));
fprintf('  r (annual)   = %.3f%%\n', 100*Params.r);
fprintf('  tau          = %.3f\n', Params.tau);
fprintf('  theta        = %.3f\n', Params.theta);

end % doPart(1)

% ---------------------------------------------------------------------------
% Boskin a2 grid. Built OUTSIDE the doPart gates because parts 4 (Boskin
% w/medical) and 6 (age profiles) both need n_a_B and a_grid_B even when part 2
% is not being run in this job.
% ---------------------------------------------------------------------------
% Grid for a2: same k^2.35 curvature as a1.
% a1 must match the U.S. economy's asset grid exactly: the two steady states are
% compared agent-by-agent in the CEV calculation, so the private-asset state has
% to mean the same thing in both.
n_a_B=[n_a(1),50];
%
% a2max was 50, which truncated the agents the paper is about: at the first solve
% 21.6% of the TOP ability type's own mass sat exactly on the ceiling (type 20:
% 13.6%, type 19: 5.2%, type 18: 1.1%). Aggregate mass there was only 1.6e-4
% because those types have tiny population weights -- but the CEV is computed PER
% TYPE and Fig 3b plots it AGAINST ability, so the truncation lands precisely on
% the high-ability "big winners" (p.516) whose 10-15% consumption gain is the
% headline. The paper does not state its own ceiling (the appendix gives only "an
% upper bound imposed on the asset space"), so 50 was an assumption, not a target.
%
% It then went to 200 with 75 a2 points. Both numbers are revised here:
%  - 75 -> 50 for memory. Under idioShocks==1 the EV block of the ExpAssetz
%    DC1_GI1 raw holds [N_d2*N_a1,N_a2,N_z,N_z] arrays that lowmemory does not
%    shrink: 6.71 GB each at N_a2=75, which OOMed on 2026-10-15 (survived the EV
%    line, died on the logical-index line two statements later). N_a2=50 is
%    4.47 GB, against 2.24 GB at the N_a2=25 size that solves cleanly.
%  - 200 -> 150 because 197.3 was a NO-IDIO number. a2 is hard-bounded, not just
%    improbably large: the z grid tops out at exp(4*sigma_y1), which is 7.99 here
%    versus 11.6 when sigma_y1^2=0.376. Drawing the top gridpoint every year AND
%    working l=1 every year gives a2=98.0 at the part-1 prices (r=5.83%), so
%    P(a2>100)=0 exactly. Under l==1 the 1-1e-10 quantile is 97.8 and the 99th
%    percentile only 46.9.
%    The ceiling is sensitive to r, through 45 years of compounding at
%    (1+r)/(1+g): 79.6 at r=5%, 98.0 at 5.83%, 133.0 at 7%, 174.2 at 8%. The
%    Boskin stage solves for its own r, and Table IV puts Boskin b=0 at 5.6%
%    (below the US system's 5.8%), so 150 holds to r~7.3% -- a wide margin.
%    150 also buys ~11% finer spacing near the a2 mean (2.806) than 200, since
%    spacing scales as a2max^0.4255 under this curvature. a2max=100 would be
%    finer still but binds above r~5.9%, too close to our own 5.83%.
Params.a2max=150;
k2_idx=(0:n_a_B(2)-1)';
a2_grid=Params.a2max*(k2_idx/(n_a_B(2)-1)).^2.35;
a_grid_B=[asset_grid;a2_grid];

%% ===========================================================================
%% Phase 2: Boskin proposal, one GE per floor level (Table III/IV rows 2-5)
%% ===========================================================================
% State swap: ebar -> a2 (accumulated tax value). theta fixed at US value;
% proportionality factor C is endogenous and ensures PAYG SS balance.

if doPart(2)==1
fprintf('\n=== HV1999 Boskin proposal, floor levels %s [%s] ===\n', ...
    mat2str(Params.bfloorfracs),specstr);

% HV1999 p.514: "The social security tax rate theta and the path of government
% consumption are set equal to those values obtained under the U.S. system."
% doPart(1) sets Params.theta from its own p_eqm, but this stage is meant to be
% runnable as its own job, so recover it from disk when doPart(1) did not just
% run. Without this a doPart=[0,1,0] job would silently hold theta at the
% initial guess (0.10) rather than the U.S. equilibrium value.
if doPart(1)~=1
    if ~exist(sprintf('./SavedOutput/HV1999_USsystem_%s.mat',specstr),'file')
        error('doPart(2) needs theta from the U.S. system; run doPart(1) first')
    end
    S=load(sprintf('./SavedOutput/HV1999_USsystem_%s.mat',specstr),'p_eqm');
    Params.theta=S.p_eqm.theta; clear S
end
fprintf('Boskin: theta held at the U.S. value %.4f\n',Params.theta);



% Initial guess for C: pick something positive; GE will solve it.
Params.C=0.05;

% aprimeFn, ReturnFn, newborn distribution and FnsToEvaluate for Boskin.
% Same noidio/idio split as the U.S. system: log productivity is a per-type
% parameter in one and the z state in the other, in the same argument slot.
vfoptions_B=vfoptions;
simoptions_B=simoptions;
if idioShocks==0
    vfoptions_B.experienceasset=1; simoptions_B.experienceasset=1;
    vfoptions_B.aprimeFn=@(l,a2,agej,Jr,g,r,w,theta,kappa_j,logabilityi) ...
        HuggettVentura1999_a2primeFn(l,a2,logabilityi,agej,Jr,g,r,w,theta,kappa_j);
else
    vfoptions_B.experienceassetz=1; simoptions_B.experienceassetz=1;
    vfoptions_B.aprimeFn=@(l,a2,z,agej,Jr,g,r,w,theta,kappa_j) ...
        HuggettVentura1999_a2primeFn(l,a2,z,agej,Jr,g,r,w,theta,kappa_j);
end
simoptions_B.aprimeFn=vfoptions_B.aprimeFn;
simoptions_B.a_grid=a_grid_B;
simoptions_B.d_grid=d_grid;

if idioShocks==0
    ReturnFn_B=@(l,a1prime,a1,a2,logabilityi,agej,Jr,sigma,nu,g,r,w,tau,theta,kappa_j,b_lower,C) ...
        HuggettVentura1999_Boskin_ReturnFn(l,a1prime,a1,a2,logabilityi,agej,Jr,sigma,nu,g,r,w,tau,theta,kappa_j,b_lower,C);
    jequaloneDist_B=zeros(n_a_B);
    jequaloneDist_B(1,1)=1;
    FnsToEvaluate_B.K        =@(l,a1prime,a1,a2) a1;
    FnsToEvaluate_B.L        =@(l,a1prime,a1,a2,agej,Jr,kappa_j,logabilityi) (agej<Jr)*exp(logabilityi)*kappa_j*l;
    FnsToEvaluate_B.Earnings =@(l,a1prime,a1,a2,agej,Jr,kappa_j,w,logabilityi) (agej<Jr)*exp(logabilityi)*kappa_j*l*w;
    FnsToEvaluate_B.AccBeq   =@(l,a1prime,a1,a2,sj,r,tau,g) (1-sj)*a1prime*(1+g)*(1+r*(1-tau));
    FnsToEvaluate_B.SSBenefits=@(l,a1prime,a1,a2,agej,Jr,g,b_lower,C) ...
        (agej>=Jr)*max(b_lower,C*a2/((1+g)^(agej-Jr)));
    % The annuity component alone, i.e. benefits BEFORE the floor is applied.
    % theta finances this part; the excess of SSBenefits over it is the floor
    % top-up, which the general income tax finances (see TaxRate below).
    FnsToEvaluate_B.SSAnnuity=@(l,a1prime,a1,a2,agej,Jr,g,C) ...
        (agej>=Jr)*C*a2/((1+g)^(agej-Jr));
    FnsToEvaluate_B.A2       =@(l,a1prime,a1,a2) a2;
    FnsToEvaluate_B.HoursWorked=@(l,a1prime,a1,a2,agej,Jr) (agej<Jr)*l;
else
    ReturnFn_B=@(l,a1prime,a1,a2,z,agej,Jr,sigma,nu,g,r,w,tau,theta,kappa_j,b_lower,C) ...
        HuggettVentura1999_Boskin_ReturnFn(l,a1prime,a1,a2,z,agej,Jr,sigma,nu,g,r,w,tau,theta,kappa_j,b_lower,C);
    jequaloneDist_B=zeros([n_a_B,n_z]);
    jequaloneDist_B(1,1,:)=reshape(jequaloneDistz,[1,1,n_z]);
    FnsToEvaluate_B.K        =@(l,a1prime,a1,a2,z) a1;
    FnsToEvaluate_B.L        =@(l,a1prime,a1,a2,z,agej,Jr,kappa_j) (agej<Jr)*exp(z)*kappa_j*l;
    FnsToEvaluate_B.Earnings =@(l,a1prime,a1,a2,z,agej,Jr,kappa_j,w) (agej<Jr)*exp(z)*kappa_j*l*w;
    FnsToEvaluate_B.AccBeq   =@(l,a1prime,a1,a2,z,sj,r,tau,g) (1-sj)*a1prime*(1+g)*(1+r*(1-tau));
    FnsToEvaluate_B.SSBenefits=@(l,a1prime,a1,a2,z,agej,Jr,g,b_lower,C) ...
        (agej>=Jr)*max(b_lower,C*a2/((1+g)^(agej-Jr)));
    % The annuity component alone, i.e. benefits BEFORE the floor is applied.
    % theta finances this part; the excess of SSBenefits over it is the floor
    % top-up, which the general income tax finances (see TaxRate below).
    FnsToEvaluate_B.SSAnnuity=@(l,a1prime,a1,a2,z,agej,Jr,g,C) ...
        (agej>=Jr)*C*a2/((1+g)^(agej-Jr));
    FnsToEvaluate_B.A2       =@(l,a1prime,a1,a2,z) a2;
    FnsToEvaluate_B.HoursWorked=@(l,a1prime,a1,a2,z,agej,Jr) (agej<Jr)*l;
end

% GE: tau, r, w, C, b_lower, ybar_mean endogenous; theta fixed at US value.
% b_lower is an unknown because the floor is specified as a multiple of output
% per person (b = frac*Y) and Y is endogenous -- exactly the same construction as
% b_common=0.1242*Y in the U.S. system stage.
GEPriceParamNames_B={'r','w','tau','C','b_lower','ybar_mean_unwage'};

GeneralEqmEqns_B.CapitalMarket=@(r,A,alpha,K,L,delta) r-(alpha*A*K^(alpha-1)*L^(1-alpha)-delta);
GeneralEqmEqns_B.LaborMarket  =@(w,A,alpha,K,L)       w-(1-alpha)*A*K^alpha*L^(-alpha);
% Government budget. With a positive floor the social security system is no
% longer self-financing: theta covers only the annuity, and the floor top-up
% (SSBenefits-SSAnnuity, i.e. sum of max(0, b_lower-annuity) over retirees) is a
% transfer financed out of general revenue. HV1999 footnote 23 is explicit about
% this -- "transfers needed to finance b=0.35Y are on the order of 1.3-1.4% of
% GDP in the no idiosyncratic case. This results in income tax rates going from
% 18.7% in the zero floor situation to 22.3% for b=0.35Y". Our b=0 run gave
% tau=0.1903 against their 18.7%, so the zero-floor end already lines up; the
% transfer term is what has to carry tau up to ~22% at the top floor.
% At b_lower=0 the top-up is identically zero and this reduces to the old
% equation, so the frac=0 case is a regression check on this whole restructure.
GeneralEqmEqns_B.TaxRate      =@(tau,r,K,w,L,AccBeq,GoverY,Y,SSBenefits,SSAnnuity) ...
    tau*(r*K+w*L)+AccBeq-GoverY*Y-(SSBenefits-SSAnnuity);
% The floor, as a multiple of output per person.
GeneralEqmEqns_B.bFloor       =@(b_lower,Y,bfloorfrac) b_lower-bfloorfrac*Y;
% SS balance, written as "price minus target" in C rather than as a budget gap.
% Here theta is FIXED at its U.S. value and the price being solved is C, the
% annuity rate. The budget-gap form theta*w*L-SSBenefits is DECREASING in C
% (benefits are C*sum(a2/(1+g)^(j-Jr)) over retirees), so the add=0 shooting
% step p <- p - factor*residual would push C the wrong way and diverge.
% Dividing through gives C - C*theta*w*L/SSBenefits, which has the same zero,
% is increasing in C, and is self-scaling: with benefits linear in C the target
% C*theta*w*L/SSBenefits is exactly theta*w*L/sum(a2/(1+g)^(j-Jr)), the full
% Newton step, so factor=0.5 is a clean half-step with no need to know that sum.
% NB this divides by SSAnnuity, which is positive whenever retirees hold any
% a2; with C>0 that holds, but it is the one term to check first if the Boskin
% solve misbehaves early.
%
% It is SSAnnuity and not SSBenefits that belongs here. HV1999 p.514: "The
% proportionality factor C is then set so that benefit payments to current
% retirees equal current social security tax payments". With a positive floor,
% total benefits exceed the annuity, so pinning C off SSBenefits would instead
% force TOTAL benefits down to theta*w*L -- driving C toward zero as the floor
% rises and leaving no transfer, which contradicts footnote 23. Pinning C off
% the annuity keeps theta financing the annuity on a pay-as-you-go basis (the
% paper's "financing of benefits is always on a pay-as-you-go basis", p.515) and
% leaves the floor top-up to the income tax. SSAnnuity is linear in C just as
% SSBenefits was at a zero floor, so the self-scaling full-Newton property of
% this form is unchanged.
GeneralEqmEqns_B.SSBalance    =@(C,theta,w,L,SSAnnuity) C-C*theta*w*L/SSAnnuity;
% ybar_mean_unwage stays a GE unknown under Boskin too: it no longer enters the
% benefit (which is C*a2), but it still defines the ebar units and so must be
% the same object in both economies for the CEV comparison to be like-for-like.
GeneralEqmEqns_B.YbarMean     =@(ybar_mean_unwage,Earnings,w,workshare) ybar_mean_unwage-Earnings/(w*workshare);

heteroagentoptions_B=heteroagentoptions;
heteroagentoptions_B.intermediateEqns.Y=@(A,alpha,K,L) A*K^alpha*L^(1-alpha);
% Override the inherited howtoupdate: Boskin solves 6 prices, not the U.S.
% system's 6 identical ones (theta is fixed here and C is new; b_lower plays the
% role b_common plays there, and like b_common it is just frac*Y so it takes the
% same full step).
heteroagentoptions_B.fminalgo9.howtoupdate={ ...
    'CapitalMarket', 'r',   0, 0.3; ...
    'LaborMarket',   'w',   0, 0.5; ...
    'TaxRate',       'tau', 0, 0.5; ...
    'SSBalance',     'C',   0, 0.5; ...
    'bFloor',        'b_lower', 0, 1.0; ...
    'YbarMean', 'ybar_mean_unwage', 0, 0.5};

% One GE per floor level, over whichever levels doFloors selects. Everything
% above this point is floor-independent, so it is built once. theta is untouched
% throughout (it is not in GEPriceParamNames_B), so it stays at the U.S. value
% as the paper requires.
for bb=find(doFloorsB)

Params.bfloorfrac=Params.bfloorfracs(bb);
% Warm start. When the levels run in sequence the prices carried over from the
% previous iteration are already a good guess. When a single level is run on its
% own, take the prices from the nearest lower level that is on disk instead of
% starting cold -- levels 1-3 each converged in 8-9 Anderson iterations from a
% warm start, and a cold start at a high floor is the case most likely to wander.
for pb=bb-1:-1:1
    pf=sprintf('./SavedOutput/HV1999_Boskin_bfloor%03d_%s.mat', ...
        round(1000*Params.bfloorfracs(pb)),specstr);
    if exist(pf,'file')
        SW=load(pf,'p_eqm_B');
        for pp=1:length(GEPriceParamNames_B)
            Params.(GEPriceParamNames_B{pp})=SW.p_eqm_B.(GEPriceParamNames_B{pp});
        end
        clear SW
        fprintf('Warm start for floor level %d from %s\n',bb,pf);
        break
    end
end
% b_lower itself is pinned by the GE condition to bfloorfrac*Y; this guess only
% needs the right order of magnitude. Y is around 0.70-0.74 in both specs.
if exist('Y_B','var')
    Params.b_lower=Params.bfloorfrac*Y_B;
else
    Params.b_lower=Params.bfloorfrac*0.72;
end
fprintf('\n--- Boskin floor level %d of %d: b = %.2f*Y ---\n', ...
    bb,length(Params.bfloorfracs),Params.bfloorfrac);

GEstagelabel=sprintf('Boskin bfloor=%.2fY',Params.bfloorfrac);
% Memory report immediately before the GE call. The value-function iteration
% OOMed here on 2026-10-09 under idioShocks==1, inside
% ValueFnIter_FHorz_ExpAssetz_DC1_GI1_nod1_raw, where five arrays of shape
% [N_d2*N_a1, N_a2, N_z, N_z] are live at once (Vlower, Vupper, aprimeProbs,
% skipinterp, EV). whos shows what THIS workspace is holding; gpuDevice shows
% what is actually free on the card, which is the number that matters -- a
% workspace variable that is a gpuArray occupies the card, and HV2000 was once
% OOMed by exactly that (stale V/Policy/StationaryDist left from a prior stage).
fprintf('\n---- memory before the %s GE solve ----\n',GEstagelabel);
try
    gd=gpuDevice;
    fprintf('GPU %s: %.2f GB free of %.2f GB total\n', ...
        gd.Name,gd.AvailableMemory/1e9,gd.TotalMemory/1e9);
catch ME
    fprintf('gpuDevice() unavailable: %s\n',ME.message);
end
whos
fprintf('---- end memory report ----\n\n');
[p_eqm_B,GeneralEqmCondn_B]=HeteroAgentStationaryEqm_Case1_FHorz_PType( ...
    n_d,n_a_B,n_z,N_j,Names_i,[],pi_z,d_grid,a_grid_B,z_grid,jequaloneDist_B, ...
    ReturnFn_B,FnsToEvaluate_B,GeneralEqmEqns_B,Params,DiscountFactorParamNames, ...
    AgeWeightsParamNames,PTypeDistParamNames,GEPriceParamNames_B, ...
    heteroagentoptions_B,simoptions_B,vfoptions_B);

for pp=1:length(GEPriceParamNames_B)
    Params.(GEPriceParamNames_B{pp})=p_eqm_B.(GEPriceParamNames_B{pp});
end

[V_B,Policy_B]=ValueFnIter_Case1_FHorz_PType(n_d,n_a_B,n_z,N_j,Names_i,d_grid,a_grid_B,z_grid,pi_z,ReturnFn_B,Params,DiscountFactorParamNames,vfoptions_B);
StationaryDist_B=StationaryDist_Case1_FHorz_PType(jequaloneDist_B,AgeWeightsParamNames,PTypeDistParamNames,Policy_B,n_d,n_a_B,n_z,N_j,Names_i,pi_z,Params,simoptions_B);
AggVars_B=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist_B,Policy_B,FnsToEvaluate_B,Params,n_d,n_a_B,n_z,N_j,Names_i,d_grid,a_grid_B,z_grid,simoptions_B);

%% Distributional statistics (same definitions as the U.S. system, Boskin benefit)
FnsToEvaluate_stats_B=FnsToEvaluate_B;
% AtFloor: fraction of retired agents at the floor (HV1999 Table III last
% column). With b_lower=0 the annuity is always >=0 so this is 0 by
% construction; it is here so the column is populated once b_lower>0 runs.
if idioShocks==0
    FnsToEvaluate_stats_B.Income=@(l,a1prime,a1,a2,agej,Jr,r,kappa_j,w,g,b_lower,C,logabilityi) ...
        (agej<Jr)*exp(logabilityi)*kappa_j*l*w + r*a1 + (agej>=Jr)*max(b_lower,C*a2/((1+g)^(agej-Jr)));
    FnsToEvaluate_stats_B.Wage=@(l,a1prime,a1,a2,agej,Jr,kappa_j,w,logabilityi) w*exp(logabilityi)*kappa_j;
    FnsToEvaluate_stats_B.AtFloor=@(l,a1prime,a1,a2,agej,Jr,g,b_lower,C) ...
        (agej>=Jr)*(b_lower>=C*a2/((1+g)^(agej-Jr)));
else
    FnsToEvaluate_stats_B.Income=@(l,a1prime,a1,a2,z,agej,Jr,r,kappa_j,w,g,b_lower,C) ...
        (agej<Jr)*exp(z)*kappa_j*l*w + r*a1 + (agej>=Jr)*max(b_lower,C*a2/((1+g)^(agej-Jr)));
    FnsToEvaluate_stats_B.Wage=@(l,a1prime,a1,a2,z,agej,Jr,kappa_j,w) w*exp(z)*kappa_j;
    FnsToEvaluate_stats_B.AtFloor=@(l,a1prime,a1,a2,z,agej,Jr,g,b_lower,C) ...
        (agej>=Jr)*(b_lower>=C*a2/((1+g)^(agej-Jr)));
end
simoptions_stats_B=simoptions_B;
if idioShocks==0
    simoptions_stats_B.conditionalrestrictions.Under65=@(l,a1prime,a1,a2,agej,Jr) (agej<Jr);
    simoptions_stats_B.conditionalrestrictions.Retired=@(l,a1prime,a1,a2,agej,Jr) (agej>=Jr);
else
    simoptions_stats_B.conditionalrestrictions.Under65=@(l,a1prime,a1,a2,z,agej,Jr) (agej<Jr);
    simoptions_stats_B.conditionalrestrictions.Retired=@(l,a1prime,a1,a2,z,agej,Jr) (agej>=Jr);
end
AllStats_B=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_B,Policy_B,FnsToEvaluate_stats_B,Params,n_d,n_a_B,n_z,N_j,Names_i,d_grid,a_grid_B,z_grid,simoptions_stats_B);
% Same spec split as the US system above: Table III b=0.0Y has Gini 0.39, Table
% IV b=0.0Y has 0.43. The 'retired at floor' entry is 0.0 at b=0 in both tables.
if idioShocks==0
    tableRef_B='Table III';
else
    tableRef_B='Table IV';
end
% Per-level print carries no paper reference: the floor-by-floor comparison
% against %s is printed as a block after the loop, where all four rows are in
% hand at once.
fprintf('\nDistributional checks (Boskin, b=%.2fY, b_lower=%g) [vs %s]:\n',Params.bfloorfrac,Params.b_lower,tableRef_B);
fprintf('  Income Gini (all agents)   = %.3f\n',gather(AllStats_B.Income.Gini));
fprintf('  Wage Gini (under 65)       = %.3f\n',gather(AllStats_B.Under65.Wage.Gini));
fprintf('  %% retired at floor         = %.1f\n',100*gather(AllStats_B.Retired.AtFloor.Mean));

%% CheckTopOfGrid: mass at the top of each endogenous grid.
% StationaryDist_B is [n_a_B(1) assets a1, n_a_B(2) accumulated SS taxes a2, n_z, N_j].
% Unlike the U.S. system's ebar grid, BOTH dimensions here are ordinary
% discretizations of continuous states with a chosen ceiling (amax, a2max), so
% mass at the top of either is a numerical artifact and both are warned on.
% Pool over permanent types, as for the U.S. system above.
pooled_B=zeros(n_a_B);
for ii=1:N_i
    di=gather(StationaryDist_B.(Names_i{ii}));
    pooled_B=pooled_B+StationaryDist_B.ptweights(ii)*sum(di,3:ndims(di)); % sum over age AND, under idioShocks==1, over z -- leaving the two endogenous dims
end
mass_by_a1_B=sum(pooled_B,2);
mass_by_a2_B=sum(pooled_B,1);
CheckTopOfGrid_B.MassAtTop10gridptsAssets=gather(sum(mass_by_a1_B(end-9:end)));
CheckTopOfGrid_B.MassAtTop5gridptsAssets =gather(sum(mass_by_a1_B(end-4:end)));
CheckTopOfGrid_B.MassAtTop2gridptsAssets =gather(sum(mass_by_a1_B(end-1:end)));
CheckTopOfGrid_B.MassAtTop10gridptsA2    =gather(sum(mass_by_a2_B(end-9:end)));
CheckTopOfGrid_B.MassAtTop5gridptsA2     =gather(sum(mass_by_a2_B(end-4:end)));
CheckTopOfGrid_B.MassAtTop2gridptsA2     =gather(sum(mass_by_a2_B(end-1:end)));
fprintf('\nCheckTopOfGrid (Boskin), amax=%g a2max=%g, threshold %.0e:\n',Params.amax,Params.a2max,checkTol);
fprintf('  mass in top 10 asset gridpts : %.4e\n',CheckTopOfGrid_B.MassAtTop10gridptsAssets);
fprintf('  mass in top  5 asset gridpts : %.4e\n',CheckTopOfGrid_B.MassAtTop5gridptsAssets);
fprintf('  mass in top  2 asset gridpts : %.4e\n',CheckTopOfGrid_B.MassAtTop2gridptsAssets);
fprintf('  mass in top 10 a2 gridpts    : %.4e\n',CheckTopOfGrid_B.MassAtTop10gridptsA2);
fprintf('  mass in top  5 a2 gridpts    : %.4e\n',CheckTopOfGrid_B.MassAtTop5gridptsA2);
fprintf('  mass in top  2 a2 gridpts    : %.4e\n',CheckTopOfGrid_B.MassAtTop2gridptsA2);
if CheckTopOfGrid_B.MassAtTop10gridptsAssets>checkTol
    warning('Boskin: %.4e mass in top 10 asset gridpts (>%.0e); widen Params.amax',CheckTopOfGrid_B.MassAtTop10gridptsAssets,checkTol);
end
if CheckTopOfGrid_B.MassAtTop10gridptsA2>checkTol
    warning('Boskin: %.4e mass in top 10 a2 gridpts (>%.0e); widen Params.a2max',CheckTopOfGrid_B.MassAtTop10gridptsA2,checkTol);
end

Y_B=Params.A*AggVars_B.K.Mean^Params.alpha*AggVars_B.L.Mean^(1-Params.alpha);
KoverY_B=AggVars_B.K.Mean/Y_B;

% One .mat per floor level. The tag is the floor fraction in per-mille so the
% names sort in the same order as the Table III/IV rows: bfloor000, bfloor150,
% bfloor250, bfloor350.
bftag=sprintf('bfloor%03d',round(1000*Params.bfloorfrac));

% Summary row for the Table III/IV block. Stored INSIDE this level's own .mat
% rather than accumulated in the workspace, because a job may solve only a
% subset of floor levels; doPart(7) assembles the table from the four files.
Row.bfloorfrac =Params.bfloorfrac;
Row.b_lower    =Params.b_lower;
Row.KoverY     =KoverY_B;
Row.Y          =Y_B;
Row.L          =gather(AggVars_B.L.Mean);
Row.HoursWorked=gather(AggVars_B.HoursWorked.Mean)/sum(Params.mewj(Params.agej<Params.Jr));
Row.r          =Params.r;
Row.tau        =Params.tau;
Row.C          =Params.C;
Row.IncomeGini =gather(AllStats_B.Income.Gini);
Row.AtFloor    =gather(AllStats_B.Retired.AtFloor.Mean);
Row.Transfer   =gather(AggVars_B.SSBenefits.Mean-AggVars_B.SSAnnuity.Mean);

fprintf('\nBoskin b=%.2fY [%s] headline aggregates:\n',Params.bfloorfrac,specstr);
fprintf('  b_lower      = %.4f\n', Params.b_lower);
fprintf('  K/Y          = %.3f\n', KoverY_B);
fprintf('  L            = %.3f\n', AggVars_B.L.Mean);
fprintf('  HoursWorked  = %.3f\n', Row.HoursWorked);
fprintf('  r (annual)   = %.3f%%\n', 100*Params.r);
fprintf('  tau          = %.3f\n', Params.tau);
fprintf('  theta        = %.3f (fixed at US value)\n', Params.theta);
fprintf('  C            = %.4f\n', Params.C);
fprintf('  floor transfer / Y = %.4f  (HV1999 fn23: ~0.013-0.014 at b=0.35Y, noidio)\n', ...
    Row.Transfer/Y_B);

save(sprintf('./SavedOutput/HV1999_Boskin_%s_%s.mat',bftag,specstr),'Params','V_B','Policy_B','StationaryDist_B','AggVars_B','AllStats_B','p_eqm_B','GeneralEqmCondn_B','CheckTopOfGrid_B');

% Free the large gpuArrays before the next level's solve. They are on the card
% until reassigned part-way through the next iteration (V_B 0.17 GB, Policy_B
% 0.68 GB, StationaryDist_B 0.17 GB at these grid sizes -- which matches the
% 1.3 GB drop seen between the level-1 and level-2 memory reports). Holding them
% wastes a GB and, worse, leaves three large blocks in the middle of the heap
% across every subsequent solve. Everything needed later is already in the .mat.
clear V_B Policy_B StationaryDist_B AggVars_B AllStats_B Row

end % floor level loop (bb)

% The Table III/IV comparison is assembled in doPart(7) from the per-level .mat
% files, not here: a job may solve only the floor levels doFloors selects, so
% the full table is not in hand at the end of this stage.

end % doPart(2)

%% ===========================================================================
%% Phase 3: Compensating variations by ability, all floor levels (Fig 3a/3b)
%% ===========================================================================
% CEV_i: percentage change in lifetime consumption such that an ability-i
% newborn is indifferent between U.S. system and Boskin. For CRRA over
% (c^nu(1-l)^{1-nu}), scaling c by (1+x) scales utility by (1+x)^{nu(1-sigma)}.
% So x = (V_US/V_B)^{1/(nu(1-sigma))} - 1.
%
% Newborn states: a=0, ebar=0 (US) / a1=0, a2=0 (Boskin); ability z indexed by i.
if doPart(3)==1

% Each stage is meant to be runnable as its own cluster job, so do not assume
% V_US is still in the workspace from doPart(1); read it back when it is not.
if ~exist('V_US','var')
    if ~exist(sprintf('./SavedOutput/HV1999_USsystem_%s.mat',specstr),'file')
        error('doPart(3) needs ./SavedOutput/HV1999_USsystem_%s.mat; run doPart(1) first',specstr)
    end
    S=load(sprintf('./SavedOutput/HV1999_USsystem_%s.mat',specstr),'V_US'); V_US=S.V_US; clear S
end

% Ability axis and the newborn ability distribution. Under idioShocks==0 each
% ability is a permanent type, so V is a struct and the newborn value is that
% type's V at (a=0, ebar=0, age 1); the type weights are abilityi_dist. Under
% idioShocks==1 there is a single type and ability at birth is the z index, so
% the newborn values are a slice along the z dimension of the one V and the
% weights are jequaloneDistz. Either way the index is ability at birth, which is
% what HV1999 Fig 3 plots on its horizontal axis.
if idioShocks==0
    nAbility=N_i;
    abilityAtBirth=logability_grid;
    abilityWeights=Params.abilityi_dist(:);
else
    nAbility=n_z;
    abilityAtBirth=z_grid;
    abilityWeights=jequaloneDistz(:);
end

% Newborn value under the U.S. system (common to every floor level).
V0_US=zeros(nAbility,1);
if idioShocks==0
    for ii=1:nAbility
        vUS=V_US.(Names_i{ii}); V0_US(ii)=gather(vUS(1,1,1));
    end
else
    vUS=V_US.(Names_i{1});
    % V is [n_a(1) x n_a(2) x n_z x N_j]; newborn is (a=0, ebar/a2=0, z, age 1).
    V0_US=gather(reshape(vUS(1,1,:,1),[nAbility,1]));
end

expo=Params.nu*(1-Params.sigma);
nFloors=length(Params.bfloorfracs);
CEVall=zeros(nAbility,nFloors);
V0_Ball=zeros(nAbility,nFloors);
CEVagg=zeros(1,nFloors);      % HV1999 Table V: ex-ante (behind the veil) CEV
ShareLoss=zeros(1,nFloors);   % HV1999 Table V: share of newborns worse off

for bb=1:nFloors
    bftag=sprintf('bfloor%03d',round(1000*Params.bfloorfracs(bb)));
    fn=sprintf('./SavedOutput/HV1999_Boskin_%s_%s.mat',bftag,specstr);
    if ~exist(fn,'file')
        error('doPart(3) needs %s; run doPart(2) first',fn)
    end
    S=load(fn,'V_B'); V_Bbb=S.V_B; clear S
    if idioShocks==0
        for ii=1:nAbility
            vB=V_Bbb.(Names_i{ii}); V0_Ball(ii,bb)=gather(vB(1,1,1));
        end
    else
        vB=V_Bbb.(Names_i{1});
        V0_Ball(:,bb)=gather(reshape(vB(1,1,:,1),[nAbility,1]));
    end
    clear V_Bbb
    CEVall(:,bb)=real((V0_Ball(:,bb)./V0_US).^(1/expo)-1); % >0 means Boskin preferred
    % Ex-ante welfare: the CEV of expected newborn utility, i.e. the welfare
    % measure for an agent behind the veil of ignorance about ability.
    CEVagg(bb)=real((sum(abilityWeights.*V0_Ball(:,bb))/sum(abilityWeights.*V0_US))^(1/expo)-1);
    ShareLoss(bb)=sum(abilityWeights(CEVall(:,bb)<0));
end

save(sprintf('./SavedOutput/HV1999_CEV_%s.mat',specstr), ...
    'V0_US','V0_Ball','CEVall','CEVagg','ShareLoss','abilityAtBirth','abilityWeights','Names_i');

fprintf('\n===== Compensating variations, Boskin vs US [%s] (Fig 3%s) =====\n', ...
    specstr,char('b'-(idioShocks==1)));  % idio -> 'a', noidio -> 'b'
fprintf('  ability |    z     |');
for bb=1:nFloors, fprintf('  b=%.2fY',Params.bfloorfracs(bb)); end
fprintf('\n');
for ii=1:nAbility
    fprintf('    %2d    | %+6.3f  |',ii,abilityAtBirth(ii));
    fprintf('  %+7.2f',100*CEVall(ii,:));
    fprintf('\n');
end
% HV1999 Table V: the aggregate measure is "never positive for the floor levels
% considered" and "the majority of agents experience a welfare loss".
% HV1999 Table V: aggregate welfare gains and agents with welfare losses.
if idioShocks==0
    refWelf=[-3.5,-1.9,-0.9,-1.5]; refLoss=[57.9,72.5,80.5,80.5];
else
    refWelf=[-2.4,-1.3,-0.3,-0.9]; refLoss=[72.5,72.5,72.1,88.3];
end
fprintf('\n  HV1999 Table V [%s]:\n',specstr);
fprintf('    floor     ex-ante CEV (%%)      share worse off (%%)\n');
fprintf('               ours    paper         ours    paper\n');
for bb=1:nFloors
    fprintf('    b=%.2fY  %+7.3f  %+6.1f        %6.1f   %6.1f\n', ...
        Params.bfloorfracs(bb),100*CEVagg(bb),refWelf(bb),100*ShareLoss(bb),refLoss(bb));
end
if any(CEVagg>0)
    warning('HV1999 Table V says the ex-ante CEV is never positive, but floor level(s) %s give a gain', ...
        mat2str(Params.bfloorfracs(CEVagg>0)));
end

end % doPart(3)

%% ===========================================================================
%% Phase 4: Boskin proposal RETAINING the hospital/medical transfer
%% ===========================================================================
% HV1999 Fig 4 and Table VI, specified in footnote 29 (p.519). The motivation
% (p.519): under the plain Boskin proposal "median-ability agents lose the
% common hospital and medical transfer built into the U.S. system, as well as
% some of the redistribution from the concave benefit formula". This stage asks
% what happens when only the OLD-AGE component is reformed and the medical
% transfer is kept in both systems.
%
% Footnote 29: "the payroll tax in the U.S. system is split into two parts
% theta1+theta2=theta. The first part is obtained as the one that finances the
% common transfer b in the U.S. system. Once the taxes are found, they are kept
% constant in the calculations for the Boskin proposal. Notice that since theta1
% finances the common transfer, theta2 is the rate at which labor earnings
% accumulate in the social security scheme. Benefits at retirement are then
% equal to max{b_bar, b + b(a2,j)}."
if doPart(4)==1

% theta, b_common, w and L all come from the U.S. equilibrium. theta1 is the
% payroll rate that financed the common transfer there: theta1*w*L equals total
% common-transfer payments, which are b_common times the mass of retirees.
if ~exist(sprintf('./SavedOutput/HV1999_USsystem_%s.mat',specstr),'file')
    error('doPart(4) needs ./SavedOutput/HV1999_USsystem_%s.mat; run doPart(1) first',specstr)
end
S=load(sprintf('./SavedOutput/HV1999_USsystem_%s.mat',specstr),'p_eqm','AggVars_US');
Params.theta=S.p_eqm.theta;
retiredmass=sum(Params.mewj(Params.agej>=Params.Jr));
Params.theta1=S.p_eqm.b_common*retiredmass/(S.p_eqm.w*gather(S.AggVars_US.L.Mean));
Params.theta2=Params.theta-Params.theta1;
clear S
fprintf('\n=== HV1999 Boskin w/medical, floor levels %s [%s] ===\n', ...
    mat2str(Params.bfloorfracs),specstr);
fprintf('Payroll split at U.S. values: theta=%.4f = theta1 %.4f (medical) + theta2 %.4f (annuity)\n', ...
    Params.theta,Params.theta1,Params.theta2);
if Params.theta2<=0
    error('theta2=%.4f is not positive: the medical transfer alone exhausts the payroll tax, so there is nothing left to accumulate a2 with',Params.theta2)
end

% Same a2 grid and solver options as the plain Boskin stage.
n_a_M=n_a_B;
a_grid_M=a_grid_B;
Params.C=0.05;

vfoptions_M=vfoptions;
simoptions_M=simoptions;
if idioShocks==0
    vfoptions_M.experienceasset=1; simoptions_M.experienceasset=1;
    % a2 accumulates at theta2, NOT theta: footnote 29. The a2primeFn is reused
    % unchanged -- the driver simply passes theta2 into its theta slot.
    vfoptions_M.aprimeFn=@(l,a2,agej,Jr,g,r,w,theta2,kappa_j,logabilityi) ...
        HuggettVentura1999_a2primeFn(l,a2,logabilityi,agej,Jr,g,r,w,theta2,kappa_j);
else
    vfoptions_M.experienceassetz=1; simoptions_M.experienceassetz=1;
    vfoptions_M.aprimeFn=@(l,a2,z,agej,Jr,g,r,w,theta2,kappa_j) ...
        HuggettVentura1999_a2primeFn(l,a2,z,agej,Jr,g,r,w,theta2,kappa_j);
    simoptions_M.z_grid=z_grid;
end
simoptions_M.aprimeFn=vfoptions_M.aprimeFn;
simoptions_M.a_grid=a_grid_M;
simoptions_M.d_grid=d_grid;

if idioShocks==0
    ReturnFn_M=@(l,a1prime,a1,a2,logabilityi,agej,Jr,sigma,nu,g,r,w,tau,theta,kappa_j,b_lower,b_common,C) ...
        HuggettVentura1999_BoskinMed_ReturnFn(l,a1prime,a1,a2,logabilityi,agej,Jr,sigma,nu,g,r,w,tau,theta,kappa_j,b_lower,b_common,C);
    jequaloneDist_M=zeros(n_a_M);
    jequaloneDist_M(1,1)=1;
    FnsToEvaluate_M.K        =@(l,a1prime,a1,a2) a1;
    FnsToEvaluate_M.L        =@(l,a1prime,a1,a2,agej,Jr,kappa_j,logabilityi) (agej<Jr)*exp(logabilityi)*kappa_j*l;
    FnsToEvaluate_M.Earnings =@(l,a1prime,a1,a2,agej,Jr,kappa_j,w,logabilityi) (agej<Jr)*exp(logabilityi)*kappa_j*l*w;
    FnsToEvaluate_M.AccBeq   =@(l,a1prime,a1,a2,sj,r,tau,g) (1-sj)*a1prime*(1+g)*(1+r*(1-tau));
    FnsToEvaluate_M.SSBenefits=@(l,a1prime,a1,a2,agej,Jr,g,b_lower,b_common,C) ...
        (agej>=Jr)*max(b_lower,b_common+C*a2/((1+g)^(agej-Jr)));
    FnsToEvaluate_M.SSAnnuity=@(l,a1prime,a1,a2,agej,Jr,g,C) ...
        (agej>=Jr)*C*a2/((1+g)^(agej-Jr));
    FnsToEvaluate_M.A2       =@(l,a1prime,a1,a2) a2;
    FnsToEvaluate_M.HoursWorked=@(l,a1prime,a1,a2,agej,Jr) (agej<Jr)*l;
else
    ReturnFn_M=@(l,a1prime,a1,a2,z,agej,Jr,sigma,nu,g,r,w,tau,theta,kappa_j,b_lower,b_common,C) ...
        HuggettVentura1999_BoskinMed_ReturnFn(l,a1prime,a1,a2,z,agej,Jr,sigma,nu,g,r,w,tau,theta,kappa_j,b_lower,b_common,C);
    jequaloneDist_M=zeros([n_a_M,n_z]);
    jequaloneDist_M(1,1,:)=reshape(jequaloneDistz,[1,1,n_z]);
    FnsToEvaluate_M.K        =@(l,a1prime,a1,a2,z) a1;
    FnsToEvaluate_M.L        =@(l,a1prime,a1,a2,z,agej,Jr,kappa_j) (agej<Jr)*exp(z)*kappa_j*l;
    FnsToEvaluate_M.Earnings =@(l,a1prime,a1,a2,z,agej,Jr,kappa_j,w) (agej<Jr)*exp(z)*kappa_j*l*w;
    FnsToEvaluate_M.AccBeq   =@(l,a1prime,a1,a2,z,sj,r,tau,g) (1-sj)*a1prime*(1+g)*(1+r*(1-tau));
    FnsToEvaluate_M.SSBenefits=@(l,a1prime,a1,a2,z,agej,Jr,g,b_lower,b_common,C) ...
        (agej>=Jr)*max(b_lower,b_common+C*a2/((1+g)^(agej-Jr)));
    FnsToEvaluate_M.SSAnnuity=@(l,a1prime,a1,a2,z,agej,Jr,g,C) ...
        (agej>=Jr)*C*a2/((1+g)^(agej-Jr));
    FnsToEvaluate_M.A2       =@(l,a1prime,a1,a2,z) a2;
    FnsToEvaluate_M.HoursWorked=@(l,a1prime,a1,a2,z,agej,Jr) (agej<Jr)*l;
end

% GE. Seven unknowns: the plain Boskin five plus b_common (kept at 0.1242*Y in
% this economy too, since the medical transfer is retained) and b_lower.
GEPriceParamNames_M={'r','w','tau','C','b_common','b_lower','ybar_mean_unwage'};
GeneralEqmEqns_M.CapitalMarket=@(r,A,alpha,K,L,delta) r-(alpha*A*K^(alpha-1)*L^(1-alpha)-delta);
GeneralEqmEqns_M.LaborMarket  =@(w,A,alpha,K,L)       w-(1-alpha)*A*K^alpha*L^(-alpha);
% Government budget. The whole payroll tax theta*w*L is available to the social
% security system; anything it pays beyond that is a transfer out of general
% revenue, exactly as in the plain Boskin stage. Here the gap can be nonzero
% even at b_lower=0, because theta1 is frozen at its U.S. value while b_common
% is re-pinned to 0.1242*Y in the new economy's own Y.
GeneralEqmEqns_M.TaxRate      =@(tau,r,K,w,L,AccBeq,GoverY,Y,SSBenefits,theta) ...
    tau*(r*K+w*L)+AccBeq-GoverY*Y-(SSBenefits-theta*w*L);
% C is pinned so the annuity is financed by theta2 alone (footnote 29: theta2 is
% "the rate at which labor earnings accumulate in the social security scheme").
% Same price-minus-target form and the same self-scaling property as the plain
% Boskin SSBalance; see the long note on that equation above.
GeneralEqmEqns_M.SSBalance    =@(C,theta2,w,L,SSAnnuity) C-C*theta2*w*L/SSAnnuity;
GeneralEqmEqns_M.bCommon      =@(b_common,Y) b_common-0.1242*Y;
GeneralEqmEqns_M.bFloor       =@(b_lower,Y,bfloorfrac) b_lower-bfloorfrac*Y;
GeneralEqmEqns_M.YbarMean     =@(ybar_mean_unwage,Earnings,w,workshare) ybar_mean_unwage-Earnings/(w*workshare);

heteroagentoptions_M=heteroagentoptions;
heteroagentoptions_M.intermediateEqns.Y=@(A,alpha,K,L) A*K^alpha*L^(1-alpha);
heteroagentoptions_M.fminalgo9.howtoupdate={ ...
    'CapitalMarket', 'r',   0, 0.3; ...
    'LaborMarket',   'w',   0, 0.5; ...
    'TaxRate',       'tau', 0, 0.5; ...
    'SSBalance',     'C',   0, 0.5; ...
    'bCommon',       'b_common', 0, 1.0; ...
    'bFloor',        'b_lower', 0, 1.0; ...
    'YbarMean', 'ybar_mean_unwage', 0, 0.5};

for bb=find(doFloorsM)

Params.bfloorfrac=Params.bfloorfracs(bb);
% Warm start from the nearest lower floor level on disk, as in doPart(2).
for pb=bb-1:-1:1
    pf=sprintf('./SavedOutput/HV1999_BoskinMed_bfloor%03d_%s.mat', ...
        round(1000*Params.bfloorfracs(pb)),specstr);
    if exist(pf,'file')
        SW=load(pf,'p_eqm_M');
        for pp=1:length(GEPriceParamNames_M)
            Params.(GEPriceParamNames_M{pp})=SW.p_eqm_M.(GEPriceParamNames_M{pp});
        end
        clear SW
        fprintf('Warm start for floor level %d from %s\n',bb,pf);
        break
    end
end
if exist('Y_M','var')
    Params.b_lower=Params.bfloorfrac*Y_M;
    Params.b_common=0.1242*Y_M;
else
    Params.b_lower=Params.bfloorfrac*0.72;
    Params.b_common=0.1242*0.72;
end
fprintf('\n--- Boskin w/medical, floor level %d of %d: b = %.2f*Y ---\n', ...
    bb,length(Params.bfloorfracs),Params.bfloorfrac);

GEstagelabel=sprintf('BoskinMed bfloor=%.2fY',Params.bfloorfrac);
fprintf('\n---- memory before the %s GE solve ----\n',GEstagelabel);
try
    gd=gpuDevice;
    fprintf('GPU %s: %.2f GB free of %.2f GB total\n', ...
        gd.Name,gd.AvailableMemory/1e9,gd.TotalMemory/1e9);
catch ME
    fprintf('gpuDevice() unavailable: %s\n',ME.message);
end
fprintf('---- end memory report ----\n\n');

[p_eqm_M,GeneralEqmCondn_M]=HeteroAgentStationaryEqm_Case1_FHorz_PType( ...
    n_d,n_a_M,n_z,N_j,Names_i,[],pi_z,d_grid,a_grid_M,z_grid,jequaloneDist_M, ...
    ReturnFn_M,FnsToEvaluate_M,GeneralEqmEqns_M,Params,DiscountFactorParamNames, ...
    AgeWeightsParamNames,PTypeDistParamNames,GEPriceParamNames_M, ...
    heteroagentoptions_M,simoptions_M,vfoptions_M);

for pp=1:length(GEPriceParamNames_M)
    Params.(GEPriceParamNames_M{pp})=p_eqm_M.(GEPriceParamNames_M{pp});
end

[V_M,Policy_M]=ValueFnIter_Case1_FHorz_PType(n_d,n_a_M,n_z,N_j,Names_i,d_grid,a_grid_M,z_grid,pi_z,ReturnFn_M,Params,DiscountFactorParamNames,vfoptions_M);
StationaryDist_M=StationaryDist_Case1_FHorz_PType(jequaloneDist_M,AgeWeightsParamNames,PTypeDistParamNames,Policy_M,n_d,n_a_M,n_z,N_j,Names_i,pi_z,Params,simoptions_M);
AggVars_M=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist_M,Policy_M,FnsToEvaluate_M,Params,n_d,n_a_M,n_z,N_j,Names_i,d_grid,a_grid_M,z_grid,simoptions_M);

FnsToEvaluate_stats_M=FnsToEvaluate_M;
if idioShocks==0
    FnsToEvaluate_stats_M.Income=@(l,a1prime,a1,a2,agej,Jr,r,kappa_j,w,g,b_lower,b_common,C,logabilityi) ...
        (agej<Jr)*exp(logabilityi)*kappa_j*l*w + r*a1 + (agej>=Jr)*max(b_lower,b_common+C*a2/((1+g)^(agej-Jr)));
    FnsToEvaluate_stats_M.AtFloor=@(l,a1prime,a1,a2,agej,Jr,g,b_lower,b_common,C) ...
        (agej>=Jr)*(b_lower>=b_common+C*a2/((1+g)^(agej-Jr)));
else
    FnsToEvaluate_stats_M.Income=@(l,a1prime,a1,a2,z,agej,Jr,r,kappa_j,w,g,b_lower,b_common,C) ...
        (agej<Jr)*exp(z)*kappa_j*l*w + r*a1 + (agej>=Jr)*max(b_lower,b_common+C*a2/((1+g)^(agej-Jr)));
    FnsToEvaluate_stats_M.AtFloor=@(l,a1prime,a1,a2,z,agej,Jr,g,b_lower,b_common,C) ...
        (agej>=Jr)*(b_lower>=b_common+C*a2/((1+g)^(agej-Jr)));
end
simoptions_stats_M=simoptions_M;
if idioShocks==0
    simoptions_stats_M.conditionalrestrictions.Retired=@(l,a1prime,a1,a2,agej,Jr) (agej>=Jr);
else
    simoptions_stats_M.conditionalrestrictions.Retired=@(l,a1prime,a1,a2,z,agej,Jr) (agej>=Jr);
end
AllStats_M=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_M,Policy_M,FnsToEvaluate_stats_M,Params,n_d,n_a_M,n_z,N_j,Names_i,d_grid,a_grid_M,z_grid,simoptions_stats_M);

Y_M=Params.A*AggVars_M.K.Mean^Params.alpha*AggVars_M.L.Mean^(1-Params.alpha);
KoverY_M=AggVars_M.K.Mean/Y_M;

bftag=sprintf('bfloor%03d',round(1000*Params.bfloorfrac));
save(sprintf('./SavedOutput/HV1999_BoskinMed_%s_%s.mat',bftag,specstr),'Params','V_M','Policy_M','StationaryDist_M','AggVars_M','AllStats_M','p_eqm_M','GeneralEqmCondn_M');

RowM.bfloorfrac =Params.bfloorfrac;
RowM.b_lower    =Params.b_lower;
RowM.b_common   =Params.b_common;
RowM.KoverY     =KoverY_M;
RowM.Y          =Y_M;
RowM.L          =gather(AggVars_M.L.Mean);
RowM.HoursWorked=gather(AggVars_M.HoursWorked.Mean)/sum(Params.mewj(Params.agej<Params.Jr));
RowM.r          =Params.r;
RowM.tau        =Params.tau;
RowM.C          =Params.C;
RowM.IncomeGini =gather(AllStats_M.Income.Gini);
RowM.AtFloor    =gather(AllStats_M.Retired.AtFloor.Mean);

fprintf('\nBoskin w/medical b=%.2fY [%s] headline aggregates:\n',Params.bfloorfrac,specstr);
fprintf('  b_lower      = %.4f   b_common = %.4f\n', Params.b_lower, Params.b_common);
fprintf('  K/Y          = %.3f\n', KoverY_M);
fprintf('  L            = %.3f\n', AggVars_M.L.Mean);
fprintf('  HoursWorked  = %.3f\n', RowM.HoursWorked);
fprintf('  r (annual)   = %.3f%%\n', 100*Params.r);
fprintf('  tau          = %.3f\n', Params.tau);
fprintf('  C            = %.4f\n', Params.C);
fprintf('  %% retired at floor = %.1f\n', 100*RowM.AtFloor);

clear V_M Policy_M StationaryDist_M AggVars_M AllStats_M RowM

end % floor level loop (bb), Boskin w/medical

% The Table VI comparison is assembled in doPart(7) from the per-level .mat
% files, for the same reason as the Boskin table above.

end % doPart(4)

%% ===========================================================================
%% Phase 5: Compensating variations for the medical variant (Fig 4)
%% ===========================================================================
if doPart(5)==1

if ~exist('V_US','var')
    if ~exist(sprintf('./SavedOutput/HV1999_USsystem_%s.mat',specstr),'file')
        error('doPart(5) needs ./SavedOutput/HV1999_USsystem_%s.mat; run doPart(1) first',specstr)
    end
    S=load(sprintf('./SavedOutput/HV1999_USsystem_%s.mat',specstr),'V_US'); V_US=S.V_US; clear S
end
if idioShocks==0
    nAbility=N_i; abilityAtBirth=logability_grid; abilityWeights=Params.abilityi_dist(:);
else
    nAbility=n_z; abilityAtBirth=z_grid; abilityWeights=jequaloneDistz(:);
end
V0_US=zeros(nAbility,1);
if idioShocks==0
    for ii=1:nAbility
        vUS=V_US.(Names_i{ii}); V0_US(ii)=gather(vUS(1,1,1));
    end
else
    vUS=V_US.(Names_i{1});
    V0_US=gather(reshape(vUS(1,1,:,1),[nAbility,1]));
end

expo=Params.nu*(1-Params.sigma);
nFloors=length(Params.bfloorfracs);
CEVallM=zeros(nAbility,nFloors);
V0_Mall=zeros(nAbility,nFloors);
CEVaggM=zeros(1,nFloors);
ShareLossM=zeros(1,nFloors);
for bb=1:nFloors
    bftag=sprintf('bfloor%03d',round(1000*Params.bfloorfracs(bb)));
    fn=sprintf('./SavedOutput/HV1999_BoskinMed_%s_%s.mat',bftag,specstr);
    if ~exist(fn,'file')
        error('doPart(5) needs %s; run doPart(4) first',fn)
    end
    S=load(fn,'V_M'); V_Mbb=S.V_M; clear S
    if idioShocks==0
        for ii=1:nAbility
            vM=V_Mbb.(Names_i{ii}); V0_Mall(ii,bb)=gather(vM(1,1,1));
        end
    else
        vM=V_Mbb.(Names_i{1});
        V0_Mall(:,bb)=gather(reshape(vM(1,1,:,1),[nAbility,1]));
    end
    clear V_Mbb
    CEVallM(:,bb)=real((V0_Mall(:,bb)./V0_US).^(1/expo)-1);
    CEVaggM(bb)=real((sum(abilityWeights.*V0_Mall(:,bb))/sum(abilityWeights.*V0_US))^(1/expo)-1);
    ShareLossM(bb)=sum(abilityWeights(CEVallM(:,bb)<0));
end

save(sprintf('./SavedOutput/HV1999_CEVmed_%s.mat',specstr), ...
    'V0_US','V0_Mall','CEVallM','CEVaggM','ShareLossM','abilityAtBirth','abilityWeights','Names_i');

fprintf('\n===== Compensating variations, Boskin w/medical vs US [%s] (Fig 4%s) =====\n', ...
    specstr,char('b'-(idioShocks==1)));
fprintf('  ability |    z     |');
for bb=1:nFloors, fprintf('  b=%.2fY',Params.bfloorfracs(bb)); end
fprintf('\n');
for ii=1:nAbility
    fprintf('    %2d    | %+6.3f  |',ii,abilityAtBirth(ii));
    fprintf('  %+7.2f',100*CEVallM(ii,:));
    fprintf('\n');
end
% HV1999 Table VI: the same welfare measures with medical benefits held constant.
if idioShocks==0
    refWelfM=[-1.3,-1.3,-0.3,-0.4]; refLossM=[72.5,72.6,76.1,69.0];
else
    refWelfM=[-0.80,-0.81,0.16,-0.15]; refLossM=[72.6,72.6,64.5,69.0];
end
fprintf('\n  HV1999 Table VI [%s]:\n',specstr);
fprintf('    floor     ex-ante CEV (%%)      share worse off (%%)\n');
fprintf('               ours    paper         ours    paper\n');
for bb=1:nFloors
    fprintf('    b=%.2fY  %+7.3f  %+6.2f        %6.1f   %6.1f\n', ...
        Params.bfloorfracs(bb),100*CEVaggM(bb),refWelfM(bb),100*ShareLossM(bb),refLossM(bb));
end

end % doPart(5)

%% ===========================================================================
%% Phase 6: Age profiles of hours and consumption (Fig 5, Fig 6)
%% ===========================================================================
% Both figures are captioned "No idiosyncratic uncertainty", so they are only
% produced for idioShocks==0. No GE is solved here: the profiles are computed
% from the Policy and StationaryDist already saved by parts 1, 2 and 4.
%
% Fig 5: (a) average hours, U.S. vs Boskin floor=0; (b) same at floor=0.35Y;
%        (c) hours of low-ability agents (ability level 5).
% Fig 6: (a) consumption profile U.S.; (b) Boskin zero floor; (c) Boskin
%        floor=0.35Y.
if doPart(6)==1
if idioShocks~=0
    fprintf('\ndoPart(6) skipped: HV1999 Fig 5 and Fig 6 are both "No idiosyncratic\nuncertainty", so the age profiles are only built for idioShocks==0.\n');
else

% The U.S. and Boskin economies have DIFFERENT experience-asset laws of motion
% (ebar vs a2), so each needs its own simoptions. LifeCycleProfiles reconstructs
% aprime for the experience-asset dimension from that economy's own aprimeFn, so
% reusing one set across both would evaluate the Boskin profiles under the ebar
% law and silently return wrong numbers.
% Two copies of each: groupptypesforstats=1 gives the pooled profile (Fig 5a/5b,
% Fig 6 heavy line), =0 gives per-ability profiles (Fig 5c, Fig 6 light lines).
simoptions_profUS=simoptions;                   % already carries the ebar aprimeFn
simoptions_profUS.groupptypesforstats=1;
simoptions_profUSi=simoptions_profUS;
simoptions_profUSi.groupptypesforstats=0;

simoptions_profB=simoptions;
simoptions_profB.experienceasset=1;             % noidio only here, so no z variant
simoptions_profB.aprimeFn=@(l,a2,agej,Jr,g,r,w,theta,kappa_j,logabilityi) ...
    HuggettVentura1999_a2primeFn(l,a2,logabilityi,agej,Jr,g,r,w,theta,kappa_j);
simoptions_profB.a_grid=a_grid_B;
simoptions_profB.d_grid=d_grid;
simoptions_profB.groupptypesforstats=1;
simoptions_profBi=simoptions_profB;
simoptions_profBi.groupptypesforstats=0;

Profiles=struct();

% ---- U.S. system ----
S=load(sprintf('./SavedOutput/HV1999_USsystem_%s.mat',specstr),'Params','Policy_US','StationaryDist_US');
P=S.Params;   % the U.S. equilibrium's own prices
% Consumption under the U.S. system, mirroring HuggettVentura1999_USsystem_ReturnFn
% exactly (bend points 0.20/1.24, the (1+g)^(j-Jr) detrend, theta paid only
% while working). Written out rather than shared with the ReturnFn because the
% toolkit evaluates FnsToEvaluate on (d,aprime,a,...) and the ReturnFn returns
% utility, not c.
FnsProf_US.HoursWorked=@(l,aprime,a,ebar,agej,Jr) (agej<Jr)*l;
FnsProf_US.Consumption=@(l,aprime,a,ebar,agej,Jr,g,r,w,tau,theta,kappa_j,b_common,ybar_mean_unwage,logabilityi) ...
    a*(1+r*(1-tau)) ...
    + (1-tau-(agej<Jr)*theta)*(agej<Jr)*exp(logabilityi)*kappa_j*l*w ...
    + (agej>=Jr)*(b_common+w*ybar_mean_unwage*( ...
        0.90*min(ebar,0.20) ...
      + 0.32*max(0,min(ebar,1.24)-0.20) ...
      + 0.15*max(0,ebar-1.24))/((1+g)^(agej-Jr))) ...
    - (1+g)*aprime;
Profiles.US=LifeCycleProfiles_FHorz_Case1_PType(S.StationaryDist_US,S.Policy_US,FnsProf_US,P,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_profUS);
Profiles.USbyAbility=LifeCycleProfiles_FHorz_Case1_PType(S.StationaryDist_US,S.Policy_US,FnsProf_US,P,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_profUSi);
clear S P

% ---- Boskin, at each floor level ----
FnsProf_B.HoursWorked=@(l,a1prime,a1,a2,agej,Jr) (agej<Jr)*l;
FnsProf_B.Consumption=@(l,a1prime,a1,a2,agej,Jr,g,r,w,tau,theta,kappa_j,b_lower,C,logabilityi) ...
    a1*(1+r*(1-tau)) ...
    + (1-tau-(agej<Jr)*theta)*(agej<Jr)*exp(logabilityi)*kappa_j*l*w ...
    + (agej>=Jr)*max(b_lower,C*a2/((1+g)^(agej-Jr))) ...
    - (1+g)*a1prime;
for bb=1:length(Params.bfloorfracs)
    bftag=sprintf('bfloor%03d',round(1000*Params.bfloorfracs(bb)));
    fn=sprintf('./SavedOutput/HV1999_Boskin_%s_%s.mat',bftag,specstr);
    if ~exist(fn,'file')
        fprintf('doPart(6): %s missing, skipping that floor level\n',fn); continue
    end
    S=load(fn,'Params','Policy_B','StationaryDist_B');
    P=S.Params;
    Profiles.(['B_' bftag])=LifeCycleProfiles_FHorz_Case1_PType(S.StationaryDist_B,S.Policy_B,FnsProf_B,P,n_d,n_a_B,n_z,N_j,Names_i,d_grid,a_grid_B,z_grid,simoptions_profB);
    Profiles.(['BbyAbility_' bftag])=LifeCycleProfiles_FHorz_Case1_PType(S.StationaryDist_B,S.Policy_B,FnsProf_B,P,n_d,n_a_B,n_z,N_j,Names_i,d_grid,a_grid_B,z_grid,simoptions_profBi);
    clear S P
end

save(sprintf('./SavedOutput/HV1999_Profiles_%s.mat',specstr),'Profiles');
fprintf('\nAge profiles computed for Fig 5 and Fig 6 [%s]: fields %s\n', ...
    specstr,strjoin(fieldnames(Profiles)',', '));

end % idioShocks==0
end % doPart(6)

%% ===========================================================================
%% Phase 7: Figures 1-6
%% ===========================================================================
% Everything here is drawn from the saved .mat files, so this part needs no GPU
% and can be re-run on its own after any of parts 1-6.
%
% Figs 3 and 4 each have an idiosyncratic panel (a) and a no-idiosyncratic panel
% (b), so they need BOTH specs on disk. A single run is one spec, so each panel
% is drawn only if its CEV file exists and a note is printed otherwise. Run the
% driver once with idioShocks=1 and once with idioShocks=0 to get both.
%
% Figs 5 and 6 are captioned "No idiosyncratic uncertainty" and are drawn from
% the noidio profiles file alone.
if doPart(7)==1

graphdir='./SavedOutput/Graphs';
ages=Params.agejshifter+Params.agej;   % real-life age, 20..100

%% Tables III/IV (Boskin) and VI (Boskin w/medical), assembled from the
%% per-level .mat files so that floor levels solved in different jobs combine.
if idioShocks==0
    tableRefT='Table III'; refKY=[2.96,2.94,2.89,2.81]; refL=[0.411,0.407,0.401,0.391];
    refHrs=[0.330,0.324,0.315,0.307]; refr=[6.2,6.2,6.5,6.8];
    refGini=[0.39,0.40,0.40,0.41]; refFloor=[0.0,15.6,41.0,59.2];
else
    tableRefT='Table IV';  refKY=[3.11,3.08,3.02,2.95]; refL=[0.434,0.429,0.421,0.413];
    refHrs=[0.329,0.322,0.312,0.306]; refr=[5.6,5.7,6.0,6.2];
    refGini=[0.43,0.43,0.44,0.44]; refFloor=[0.0,12.7,39.1,58.8];
end
fprintf('\n===== HV1999 %s, Boskin rows [%s] =====\n',tableRefT,specstr);
fprintf('               K/Y          L          hours         r         Gini      %%at floor   transfer/Y\n');
texrows=cell(1,length(Params.bfloorfracs));
fprintf('  floor     ours paper  ours paper  ours paper  ours paper  ours paper  ours paper     ours\n');
for bb=1:length(Params.bfloorfracs)
    fn=sprintf('./SavedOutput/HV1999_Boskin_bfloor%03d_%s.mat', ...
        round(1000*Params.bfloorfracs(bb)),specstr);
    if ~exist(fn,'file')
        fprintf('  b=%.2fY   -- not yet solved (%s) --\n',Params.bfloorfracs(bb),fn); continue
    end
    % Recomputed here from the variables EVERY level's .mat contains, rather
    % than read from a summary struct: levels solved by different jobs can come
    % from different versions of this driver, and a level solved before a
    % summary field was introduced would simply not have it.
    SR=load(fn,'Params','AggVars_B','AllStats_B');
    PB=SR.Params; AV=SR.AggVars_B; AS=SR.AllStats_B; clear SR
    rK=gather(AV.K.Mean); rL=gather(AV.L.Mean);
    rY=PB.A*rK^PB.alpha*rL^(1-PB.alpha);
    fprintf('  b=%.2fY  %5.2f %5.2f  %5.3f %5.3f  %5.3f %5.3f  %4.1f %4.1f  %5.2f %5.2f  %5.1f %5.1f    %6.4f\n', ...
        PB.bfloorfrac,rK/rY,refKY(bb),rL,refL(bb), ...
        gather(AV.HoursWorked.Mean)/sum(PB.mewj(PB.agej<PB.Jr)),refHrs(bb), ...
        100*PB.r,refr(bb),gather(AS.Income.Gini),refGini(bb), ...
        100*gather(AS.Retired.AtFloor.Mean),refFloor(bb), ...
        gather(AV.SSBenefits.Mean-AV.SSAnnuity.Mean)/rY);
    texrows{bb}=sprintf('$\\underline{b}=%.2fY$ & %.2f & %.3f & %.3f & %.1f & %.2f & %.1f \\\\\n', ...
        PB.bfloorfrac,rK/rY,rL,gather(AV.HoursWorked.Mean)/sum(PB.mewj(PB.agej<PB.Jr)), ...
        100*PB.r,gather(AS.Income.Gini),100*gather(AS.Retired.AtFloor.Mean));
end

% LaTeX fragment for Table III (noidio) / Table IV (idio). Replication numbers
% ONLY: the paper's own table is shown directly above it in the write-up, so
% repeating the paper's values here would duplicate what is already on the page.
SU=load(sprintf('./SavedOutput/HV1999_USsystem_%s.mat',specstr),'Params','AggVars_US','AllStats_US');
PU=SU.Params; AU=SU.AggVars_US; ASU=SU.AllStats_US; clear SU
uK=gather(AU.K.Mean); uL=gather(AU.L.Mean); uY=PU.A*uK^PU.alpha*uL^(1-PU.alpha);
if idioShocks==0, tnum=3; else, tnum=4; end
fid=fopen(sprintf('./SavedOutput/LatexInputs/HV1999_Table%d.tex',tnum),'w');
fprintf(fid,'\\begin{tabular}{lcccccc}\n\\hline\n');
fprintf(fid,'            &       &       & Fraction  &       & Income & Percentage of \\\\\n');
fprintf(fid,'            &       &       & of time   &       & Gini   & retired agents \\\\\n');
fprintf(fid,'            & $K/Y$ & $L$   & working   & $r$   &        & at floor level \\\\\n\\hline\n');
fprintf(fid,'U.S. system & %.2f & %.3f & %.3f & %.1f & %.2f & --- \\\\\n', ...
    uK/uY,uL,gather(AU.HoursWorked.Mean)/sum(PU.mewj(PU.agej<PU.Jr)), ...
    100*PU.r,gather(ASU.Income.Gini));
fprintf(fid,'Boskin      &       &       &           &       &        &      \\\\\n');
fprintf(fid,'\\tb proposal&       &       &           &       &        &      \\\\\n');
for bb=1:length(Params.bfloorfracs)
    if ~isempty(texrows{bb}), fprintf(fid,'%s',texrows{bb}); end
end
fprintf(fid,'\\hline\n\\end{tabular}\n');
fclose(fid);
fprintf('  wrote SavedOutput/LatexInputs/HV1999_Table%d.tex\n',tnum);
fprintf('  (HV1999 fn23: the floor transfer is ~0.013-0.014 of GDP at b=0.35Y, noidio,\n');
fprintf('   with the income tax going 18.7%% -> 22.3%% across the floor range)\n');

% LaTeX fragments for Tables V and VI (the two welfare tables). Each has a
% no-idiosyncratic and an idiosyncratic column pair, so both specs must be on
% disk; whichever spec this job ran, the fragment is written only once both CEV
% files exist. Replication numbers only, as for Tables III/IV.
cevfiles={'./SavedOutput/HV1999_CEV_noidio.mat','./SavedOutput/HV1999_CEV_idio.mat'; ...
          './SavedOutput/HV1999_CEVmed_noidio.mat','./SavedOutput/HV1999_CEVmed_idio.mat'};
for tt=1:2      % tt=1 -> Table V (plain Boskin), tt=2 -> Table VI (w/medical)
    if ~exist(cevfiles{tt,1},'file') || ~exist(cevfiles{tt,2},'file')
        fprintf('  Table %d fragment not written: needs both specs (%s, %s)\n', ...
            tt+4,cevfiles{tt,1},cevfiles{tt,2}); continue
    end
    Sn=load(cevfiles{tt,1}); Si=load(cevfiles{tt,2});
    if tt==1
        wn=100*Sn.CEVagg;  ln=100*Sn.ShareLoss;  wi=100*Si.CEVagg;  li=100*Si.ShareLoss;
    else
        wn=100*Sn.CEVaggM; ln=100*Sn.ShareLossM; wi=100*Si.CEVaggM; li=100*Si.ShareLossM;
    end
    clear Sn Si
    fid=fopen(sprintf('./SavedOutput/LatexInputs/HV1999_Table%d.tex',tt+4),'w');
    fprintf(fid,'\\begin{tabular}{lcccc}\n\\hline\n');
    fprintf(fid,' & \\multicolumn{2}{c}{No idiosyncratic shocks} & \\multicolumn{2}{c}{Idiosyncratic shocks} \\\\\n');
    fprintf(fid,'\\cline{2-3} \\cline{4-5}\n');
    fprintf(fid,'Model economy & Welfare & Agents w/ & Welfare & Agents w/ \\\\\n');
    fprintf(fid,'              & gains \\%% & losses \\%% & gains \\%% & losses \\%% \\\\\n\\hline\n');
    for bb=1:length(Params.bfloorfracs)
        fprintf(fid,'$\\underline{b}=%.2fY$ & $%+.2f$ & %.1f & $%+.2f$ & %.1f \\\\\n', ...
            Params.bfloorfracs(bb),wn(bb),ln(bb),wi(bb),li(bb));
    end
    fprintf(fid,'\\hline\n\\end{tabular}\n');
    fclose(fid);
    fprintf('  wrote SavedOutput/LatexInputs/HV1999_Table%d.tex\n',tt+4);
end

% NB the paper has NO aggregates table for the medical variant -- its Table VI
% is the WELFARE table, printed by doPart(5). This block is ours, for
% comparability with Table IV.
fprintf('\n===== Boskin w/medical, aggregates [%s] (no counterpart table in HV1999) =====\n',specstr);
fprintf('  floor      K/Y      L     hours      r      Gini   %%at floor\n');
for bb=1:length(Params.bfloorfracs)
    fn=sprintf('./SavedOutput/HV1999_BoskinMed_bfloor%03d_%s.mat', ...
        round(1000*Params.bfloorfracs(bb)),specstr);
    if ~exist(fn,'file')
        fprintf('  b=%.2fY   -- not yet solved --\n',Params.bfloorfracs(bb)); continue
    end
    SR=load(fn,'Params','AggVars_M','AllStats_M');
    PM=SR.Params; AV=SR.AggVars_M; AS=SR.AllStats_M; clear SR
    rK=gather(AV.K.Mean); rL=gather(AV.L.Mean);
    rY=PM.A*rK^PM.alpha*rL^(1-PM.alpha);
    fprintf('  b=%.2fY   %5.2f  %5.3f  %5.3f  %5.2f  %5.2f  %5.1f\n', ...
        PM.bfloorfrac,rK/rY,rL, ...
        gather(AV.HoursWorked.Mean)/sum(PM.mewj(PM.agej<PM.Jr)), ...
        100*PM.r,gather(AS.Income.Gini),100*gather(AS.Retired.AtFloor.Mean));
end

% ---------------------------------------------------------------- Fig 1
% Wage profile. Params.kappa_j is the interpolated Hansen (1993) profile; the
% markers are his tabulated values at the age-group centres.
figure(1); clf
% Size and font are set HERE, before anything is drawn, not just before the
% export. A legend measures its box from the text metrics at creation time,
% so bumping the font afterwards leaves the box sized for the old font and
% the labels spill out of it -- which is exactly what happened when this was
% done just before exportgraphics: 'b=0.00Y' was clipped to 'b='.
% Setting figure defaults makes every axes, label and legend created below
% use the final size, so each legend sizes its own box correctly.
set(gcf,'Position',[100 100 820 520]);
set(gcf,'DefaultAxesFontSize',12,'DefaultTextFontSize',12,'DefaultLegendFontSize',12);
plot(ages,Params.kappa_j,'b-','LineWidth',1.5); hold on
if exist('HV1999_Hansen1993_centre','var') && exist('HV1999_Hansen1993_male','var')
    plot(HV1999_Hansen1993_centre,HV1999_Hansen1993_male,'ro', ...
        'MarkerFaceColor','r','MarkerSize',6);
    legend('Model kappa_j (interpolated)','Hansen (1993) age-group centres','Location','southwest');
else
    legend('Model kappa_j','Location','southwest');
end
xlabel('Age'); ylabel('Labor endowment in efficiency units');
title('HV1999 Fig 1: Wage profile (source: Hansen 1993)');
xlim([20 80]); grid on
% Pin the canvas size and the font before exporting. Legend boxes are sized
% from the text metrics at creation time; if the figure is then rendered at a
% different size the box no longer fits its own labels, which is why the
% legends in the first version of these figures had text spilling out of the
% box and over the curves. drawnow forces the re-measure before the capture.
drawnow;
lgs_=findobj(gcf,'Type','legend');
for q_=1:numel(lgs_)
    set(lgs_(q_),'Units','normalized');
    pq_=lgs_(q_).Position;
    lgs_(q_).Position=[pq_(1) pq_(2) pq_(3)*1.9 pq_(4)];
end
drawnow;
exportgraphics(gcf,fullfile(graphdir,'HV1999_Fig1_WageProfile.png'),'Resolution',200);

% ---------------------------------------------------------------- Fig 2
% Social security benefits and earnings-related components, as a function of
% average past earnings ebar, at the U.S. equilibrium prices. Mirrors the PIA
% schedule in HuggettVentura1999_USsystem_ReturnFn: 90% to the first bend point
% at 0.20, 32% to the second at 1.24, 15% above, with ebar capped at 2.47.
if exist(sprintf('./SavedOutput/HV1999_USsystem_%s.mat',specstr),'file')
    S=load(sprintf('./SavedOutput/HV1999_USsystem_%s.mat',specstr),'Params');
    P=S.Params; clear S
    ebar_plot=linspace(0,2.47,400)';
    PIA_plot=0.90*min(ebar_plot,0.20) ...
            +0.32*max(0,min(ebar_plot,1.24)-0.20) ...
            +0.15*max(0,ebar_plot-1.24);
    earnrel=P.w*P.ybar_mean_unwage*PIA_plot;   % at age Jr, before the (1+g) detrend
    figure(2); clf
    % Size and font are set HERE, before anything is drawn, not just before the
    % export. A legend measures its box from the text metrics at creation time,
    % so bumping the font afterwards leaves the box sized for the old font and
    % the labels spill out of it -- which is exactly what happened when this was
    % done just before exportgraphics: 'b=0.00Y' was clipped to 'b='.
    % Setting figure defaults makes every axes, label and legend created below
    % use the final size, so each legend sizes its own box correctly.
    set(gcf,'Position',[100 100 860 540]);
    set(gcf,'DefaultAxesFontSize',12,'DefaultTextFontSize',12,'DefaultLegendFontSize',12);
    plot(ebar_plot,earnrel+P.b_common,'b-','LineWidth',1.5); hold on
    plot(ebar_plot,earnrel,'r--','LineWidth',1.5);
    plot(ebar_plot,P.b_common*ones(size(ebar_plot)),'k:','LineWidth',1.2);
    xline(0.20,'Color',[.6 .6 .6]); xline(1.24,'Color',[.6 .6 .6]);
    legend('Total benefit b(x,j)','Earnings-related component','Common (medical) transfer', ...
        'Location','northwest');
    % Plain text: MATLAB's default 'tex' interpreter has no \bar, so '\bar{e}'
    % here threw "must have valid interpreter syntax" and left the label blank.
    xlabel('Average past earnings, multiples of average earnings');
    ylabel('Benefit at retirement');
    title(sprintf('HV1999 Fig 2: Social security benefits [%s]',specstr));
    grid on
    % see the note above on legend-box sizing
    drawnow;
    % Force each legend box wider than MATLAB thinks it needs. On the cluster
    % node there is no display, and text-extent measurement there comes back far
    % too small, so every legend box is sized as if its labels were nearly empty
    % and the text spills out past the border ("b=0.00Y" clipped to "b="). Four
    % attempts to make the measurement correct all failed -- enlarging the
    % canvas, a post-hoc font change, setting the font before creation, and
    % figure-level font defaults -- because the measurement itself is the thing
    % that is broken, not the scaling. Widening the box by a fixed factor is
    % empirical, but it is deterministic and it is checkable in the PNG.
    lgs_=findobj(gcf,'Type','legend');
    for q_=1:numel(lgs_)
        set(lgs_(q_),'Units','normalized');
        pq_=lgs_(q_).Position;
        lgs_(q_).Position=[pq_(1) pq_(2) pq_(3)*1.9 pq_(4)];
    end
    drawnow;
    exportgraphics(gcf,fullfile(graphdir,sprintf('HV1999_Fig2_Benefits_%s.png',specstr)),'Resolution',200);
else
    fprintf('Fig 2 skipped: no U.S. system .mat for spec %s\n',specstr);
end

% ------------------------------------------------------------- Figs 3 and 4
% Panel (a) is the idiosyncratic spec, panel (b) the no-idiosyncratic one.
specpanels={'idio','a';'noidio','b'};
for whichfig=3:4
    if whichfig==3
        fstem='CEV';    fvar='CEVall';  figname='Fig3'; figdesc='U.S. system vs. Boskin proposal';
    else
        fstem='CEVmed'; fvar='CEVallM'; figname='Fig4'; figdesc='U.S. system vs. Boskin proposal (w/medical benefits)';
    end
    figure(whichfig); clf
    % Size and font are set HERE, before anything is drawn, not just before the
    % export. A legend measures its box from the text metrics at creation time,
    % so bumping the font afterwards leaves the box sized for the old font and
    % the labels spill out of it -- which is exactly what happened when this was
    % done just before exportgraphics: 'b=0.00Y' was clipped to 'b='.
    % Setting figure defaults makes every axes, label and legend created below
    % use the final size, so each legend sizes its own box correctly.
    set(gcf,'Position',[100 100 1320 540]);
    set(gcf,'DefaultAxesFontSize',11,'DefaultTextFontSize',11,'DefaultLegendFontSize',11);
    nplotted=0;
    for sp=1:2
        thisspec=specpanels{sp,1};
        fn=sprintf('./SavedOutput/HV1999_%s_%s.mat',fstem,thisspec);
        subplot(1,2,sp)
        if ~exist(fn,'file')
            title(sprintf('(%s) %s: not yet run',specpanels{sp,2},thisspec));
            axis off
            fprintf('%s panel (%s) skipped: %s missing\n',figname,specpanels{sp,2},fn);
            continue
        end
        S=load(fn);
        CEVmat=S.(fvar);
        nAb=size(CEVmat,1);
        hold on
        for bb=1:size(CEVmat,2)
            plot(1:nAb,100*CEVmat(:,bb),'-o','LineWidth',1.3,'MarkerSize',3);
        end
        yline(0,'k:');
        xlabel('Ability level at birth'); ylabel('CEV (%)');
        if sp==1
            title(sprintf('(%s) Idiosyncratic uncertainty',specpanels{sp,2}));
        else
            title(sprintf('(%s) No idiosyncratic uncertainty',specpanels{sp,2}));
        end
        legend(arrayfun(@(f) sprintf('b=%.2fY',f),Params.bfloorfracs,'UniformOutput',false), ...
            'Location','northwest');
        xlim([1 nAb]); grid on
        nplotted=nplotted+1;
        clear S CEVmat
    end
    sgtitle(sprintf('HV1999 %s: Compensating variations. %s',figname,figdesc));
    if nplotted>0
        % see the note above on legend-box sizing
    % Wider canvas and a slightly smaller font than the other figures: the
    % ticks here are wide ("-10","-20") and at 12pt on a 1100px canvas the
    % y-label overlapped them and hid their minus signs -- "-10" rendered as
    % "10". The label is also shortened to "CEV (%)".
    % Fit each axes inside its own subplot cell using the space its decorations
    % actually need, rather than a hand-picked offset. TightInset is the margin
    % MATLAB measures for that axes' ticks, labels and title, so this adapts to
    % whatever the font and label lengths turn out to be. Two earlier attempts
    % at this by hand both failed: enlarging the canvas did nothing (subplot
    % pins the axes at a fixed normalised position, so the overlap just scales),
    % and a fixed rightward shift fixed the y-label but overwrote the band
    % sgtitle reserves, so the panels then collided with the overall title.
    % The 0.06 is that sgtitle band.
    axs_=findobj(gcf,'Type','axes');
    for kk_=1:numel(axs_)
        set(axs_(kk_),'Units','normalized');
        op_=axs_(kk_).OuterPosition; ti_=axs_(kk_).TightInset;
        axs_(kk_).Position=[op_(1)+ti_(1), op_(2)+ti_(2), ...
            max(op_(3)-ti_(1)-ti_(3),0.05), ...
            max(op_(4)-ti_(2)-ti_(4)-0.06,0.05)];
    end
    drawnow;
    % see the note at Fig 1 on forcing legend box widths
    lgs_=findobj(gcf,'Type','legend');
    for q_=1:numel(lgs_)
        set(lgs_(q_),'Units','normalized');
        pq_=lgs_(q_).Position;
        lgs_(q_).Position=[pq_(1) pq_(2) pq_(3)*1.9 pq_(4)];
    end
    drawnow;
    exportgraphics(gcf,fullfile(graphdir,sprintf('HV1999_%s_CEV.png',figname)),'Resolution',200);
    end
end

% ------------------------------------------------------------- Figs 5 and 6
% Both are no-idiosyncratic-uncertainty figures, so they read the noidio
% profiles file regardless of which spec this run was.
profn='./SavedOutput/HV1999_Profiles_noidio.mat';
if ~exist(profn,'file')
    fprintf('Figs 5 and 6 skipped: %s missing (run doPart(6) with idioShocks=0)\n',profn);
else
    S=load(profn,'Profiles'); PR=S.Profiles; clear S
    agesW=ages(Params.agej<Params.Jr);    % working ages, for the hours panels

    % ---- Fig 5: labour hours ----
    figure(5); clf
    % Size and font are set HERE, before anything is drawn, not just before the
    % export. A legend measures its box from the text metrics at creation time,
    % so bumping the font afterwards leaves the box sized for the old font and
    % the labels spill out of it -- which is exactly what happened when this was
    % done just before exportgraphics: 'b=0.00Y' was clipped to 'b='.
    % Setting figure defaults makes every axes, label and legend created below
    % use the final size, so each legend sizes its own box correctly.
    set(gcf,'Position',[100 100 1250 460]);
    set(gcf,'DefaultAxesFontSize',12,'DefaultTextFontSize',12,'DefaultLegendFontSize',12);
    % (a) floor = 0
    subplot(1,3,1); hold on
    if isfield(PR,'US') && isfield(PR,'B_bfloor000')
        hUS=gather(PR.US.HoursWorked.Mean); hB=gather(PR.B_bfloor000.HoursWorked.Mean);
        plot(agesW,hUS(Params.agej<Params.Jr),'b-','LineWidth',1.5);
        plot(agesW,hB(Params.agej<Params.Jr),'r--','LineWidth',1.5);
        legend('U.S. system','Boskin, floor = 0','Location','southwest');
    end
    xlabel('Age'); ylabel('Average labor hours'); title('(a) Floor = 0'); grid on
    % (b) hours BY ABILITY LEVEL at floor = 0.35Y. The paper's panel (b) legend
    % reads "Ability Level 5 / 10 / 15", i.e. it decomposes the profile by
    % ability rather than comparing the two systems -- panel (a) already does
    % the system comparison. An earlier version of this figure repeated the
    % system comparison here, which is not what HV1999 Fig 5(b) shows.
    subplot(1,3,2); hold on
    abset={'ability05','ability10','ability15'};
    if isfield(PR,'BbyAbility_bfloor350')
        lgb={};
        for kk=1:length(abset)
            if isfield(PR.BbyAbility_bfloor350.HoursWorked,abset{kk})
                hk=gather(PR.BbyAbility_bfloor350.HoursWorked.(abset{kk}).Mean);
                plot(agesW,hk(Params.agej<Params.Jr),'LineWidth',1.4);
                lgb{end+1}=strrep(abset{kk},'ability','Ability level ');
            end
        end
        if ~isempty(lgb), legend(lgb,'Location','southwest'); end
    end
    xlabel('Age'); ylabel('Average labor hours'); title('(b) Floor = 0.35 \times Y'); grid on
    % (c) low-ability agents, ability level 5
    subplot(1,3,3); hold on
    ab5='ability05';
    if isfield(PR,'USbyAbility') && isfield(PR.USbyAbility.HoursWorked,ab5)
        hUS5=gather(PR.USbyAbility.HoursWorked.(ab5).Mean);
        plot(agesW,hUS5(Params.agej<Params.Jr),'b-','LineWidth',1.5);
        lg={'U.S. system'};
        if isfield(PR,'BbyAbility_bfloor000') && isfield(PR.BbyAbility_bfloor000.HoursWorked,ab5)
            h05=gather(PR.BbyAbility_bfloor000.HoursWorked.(ab5).Mean);
            plot(agesW,h05(Params.agej<Params.Jr),'r--','LineWidth',1.5); lg{end+1}='Boskin, floor = 0';
        end
        if isfield(PR,'BbyAbility_bfloor350') && isfield(PR.BbyAbility_bfloor350.HoursWorked,ab5)
            h35=gather(PR.BbyAbility_bfloor350.HoursWorked.(ab5).Mean);
            plot(agesW,h35(Params.agej<Params.Jr),'g-.','LineWidth',1.5); lg{end+1}='Boskin, floor = 0.35Y';
        end
        legend(lg,'Location','southwest');
    end
    xlabel('Age'); ylabel('Average labor hours'); title('(c) Ability level 5'); grid on
    sgtitle('HV1999 Fig 5: Labor hours. No idiosyncratic uncertainty');
    % see the note above on legend-box sizing
    % Fit each axes inside its own subplot cell using the space its decorations
    % actually need, rather than a hand-picked offset. TightInset is the margin
    % MATLAB measures for that axes' ticks, labels and title, so this adapts to
    % whatever the font and label lengths turn out to be. Two earlier attempts
    % at this by hand both failed: enlarging the canvas did nothing (subplot
    % pins the axes at a fixed normalised position, so the overlap just scales),
    % and a fixed rightward shift fixed the y-label but overwrote the band
    % sgtitle reserves, so the panels then collided with the overall title.
    % The 0.06 is that sgtitle band.
    axs_=findobj(gcf,'Type','axes');
    for kk_=1:numel(axs_)
        set(axs_(kk_),'Units','normalized');
        op_=axs_(kk_).OuterPosition; ti_=axs_(kk_).TightInset;
        axs_(kk_).Position=[op_(1)+ti_(1), op_(2)+ti_(2), ...
            max(op_(3)-ti_(1)-ti_(3),0.05), ...
            max(op_(4)-ti_(2)-ti_(4)-0.06,0.05)];
    end
    drawnow;
    % see the note at Fig 1 on forcing legend box widths
    lgs_=findobj(gcf,'Type','legend');
    for q_=1:numel(lgs_)
        set(lgs_(q_),'Units','normalized');
        pq_=lgs_(q_).Position;
        lgs_(q_).Position=[pq_(1) pq_(2) pq_(3)*1.9 pq_(4)];
    end
    drawnow;
    exportgraphics(gcf,fullfile(graphdir,'HV1999_Fig5_Hours.png'),'Resolution',200);

    % ---- Fig 6: consumption ----
    % Heavy line is the population average, light lines are ability levels 5,
    % 10 and 15 -- the same three the paper's own Fig 6 legend names ("Ability
    % Level 15 / 10 / 5"). Using 5/11/21 instead put the top ability, which
    % consumes ~10x the bottom, on the same axis and flattened everything else.
    abpick={'ability05','ability10','ability15'};
    panels={'US','U.S. system'; ...
            'B_bfloor000','Boskin, zero floor'; ...
            'B_bfloor350','Boskin, floor = 0.35 \times GDP'};
    byab={'USbyAbility','BbyAbility_bfloor000','BbyAbility_bfloor350'};
    figure(6); clf
    % Size and font are set HERE, before anything is drawn, not just before the
    % export. A legend measures its box from the text metrics at creation time,
    % so bumping the font afterwards leaves the box sized for the old font and
    % the labels spill out of it -- which is exactly what happened when this was
    % done just before exportgraphics: 'b=0.00Y' was clipped to 'b='.
    % Setting figure defaults makes every axes, label and legend created below
    % use the final size, so each legend sizes its own box correctly.
    set(gcf,'Position',[100 100 1250 460]);
    set(gcf,'DefaultAxesFontSize',12,'DefaultTextFontSize',12,'DefaultLegendFontSize',12);
    for pp=1:3
        subplot(1,3,pp); hold on
        if isfield(PR,panels{pp,1})
            cAvg=gather(PR.(panels{pp,1}).Consumption.Mean);
            plot(ages,cAvg,'k-','LineWidth',2);
            lg={'All agents'};
            if isfield(PR,byab{pp})
                for kk=1:length(abpick)
                    if isfield(PR.(byab{pp}).Consumption,abpick{kk})
                        ck=gather(PR.(byab{pp}).Consumption.(abpick{kk}).Mean);
                        plot(ages,ck,'LineWidth',1.1);
                        lg{end+1}=strrep(abpick{kk},'ability','ability ');
                    end
                end
            end
            legend(lg,'Location','northwest');
        end
        xlabel('Age'); ylabel('Consumption');
        title(sprintf('(%c) %s',char('a'+pp-1),panels{pp,2})); grid on
    end
    sgtitle('HV1999 Fig 6: Consumption profiles. No idiosyncratic uncertainty');
    % see the note above on legend-box sizing
    % Fit each axes inside its own subplot cell using the space its decorations
    % actually need, rather than a hand-picked offset. TightInset is the margin
    % MATLAB measures for that axes' ticks, labels and title, so this adapts to
    % whatever the font and label lengths turn out to be. Two earlier attempts
    % at this by hand both failed: enlarging the canvas did nothing (subplot
    % pins the axes at a fixed normalised position, so the overlap just scales),
    % and a fixed rightward shift fixed the y-label but overwrote the band
    % sgtitle reserves, so the panels then collided with the overall title.
    % The 0.06 is that sgtitle band.
    axs_=findobj(gcf,'Type','axes');
    for kk_=1:numel(axs_)
        set(axs_(kk_),'Units','normalized');
        op_=axs_(kk_).OuterPosition; ti_=axs_(kk_).TightInset;
        axs_(kk_).Position=[op_(1)+ti_(1), op_(2)+ti_(2), ...
            max(op_(3)-ti_(1)-ti_(3),0.05), ...
            max(op_(4)-ti_(2)-ti_(4)-0.06,0.05)];
    end
    drawnow;
    % see the note at Fig 1 on forcing legend box widths
    lgs_=findobj(gcf,'Type','legend');
    for q_=1:numel(lgs_)
        set(lgs_(q_),'Units','normalized');
        pq_=lgs_(q_).Position;
        lgs_(q_).Position=[pq_(1) pq_(2) pq_(3)*1.9 pq_(4)];
    end
    drawnow;
    exportgraphics(gcf,fullfile(graphdir,'HV1999_Fig6_Consumption.png'),'Resolution',200);
end

fprintf('\nFigures written to %s\n',graphdir);
dd=dir(fullfile(graphdir,'*.png'));
for kk=1:length(dd), fprintf('  %s (%.0f kB)\n',dd(kk).name,dd(kk).bytes/1024); end

end % doPart(7)
%% Close out
diary off
