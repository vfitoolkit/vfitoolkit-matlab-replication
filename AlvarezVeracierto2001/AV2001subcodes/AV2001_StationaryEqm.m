function [residual,Out]=AV2001_StationaryEqm(i,Params,Grids,vfoptions,simoptions,verbose)
% Steps (1)-(10) of the AV2001 appendix for a given interest rate i.
% residual is the goods market clearing error, eq (3.14), relative to output.
% Out collects everything the tables need.
%
% Params.endogsearch=0: EXOGENOUS job-finding rate (Tables 2 and 4). Search
%   intensity is the parameter Params.eta, households have no decision variable
%   and their shocks are a plain Markov chain (n_d=0, z).
% Params.endogsearch=1: ENDOGENOUS job-finding rate (Tables 1, 3, 5). Search
%   intensity is a decision on Grids.eta_grid and the household state is a
%   semi-exogenous shock (n_d=n_eta, n_z=0, semiz) whose unemployed rows depend
%   on eta through eta^xi.
%
% Household exogenous state: the establishment states (e,s), s in {1,sbar},
% that are REACHABLE under the establishment policy (entrants start at e=0 and
% follow e'=n(e,s)), then just-fired (severance in hand), unemployed, retired,
% dead. Restricting to reachable states is exact (no household is ever employed
% elsewhere) and keeps the return matrix small: with endogenous search it is
% n_eta*n_a*n_a*N_hh doubles and the toolkit has no lowmemory path for it.

Params.i=i;
Params.r=i+Params.delta;
Params.firmdiscount=1/(1+i);

n_e=Grids.n_e; n_s=Grids.n_s; e_grid=Grids.e_grid; Q=Grids.Q;
a_grid=Grids.a_grid;
n_a=length(a_grid);
if Params.endogsearch==1
    eta_grid=Grids.eta_grid;
    n_eta=length(eta_grid);
else
    eta_grid=Params.eta; % the one exogenous value
    n_eta=1;
end

%% Step 2: establishments, and wstar from free entry
% Bracket wstar (the free entry residual is increasing in wstar), then fzero
wlo=Params.wstar_guess/1.3; whi=Params.wstar_guess*1.3;
FElo=AV2001_FreeEntryResidual(wlo,Params,Grids,vfoptions);
FEhi=AV2001_FreeEntryResidual(whi,Params,Grids,vfoptions);
bracketiter=0;
while FElo*FEhi>0 && bracketiter<30
    if FElo>0 % both positive: wstar too high
        whi=wlo; FEhi=FElo;
        wlo=wlo/1.3; FElo=AV2001_FreeEntryResidual(wlo,Params,Grids,vfoptions);
    else % both negative: wstar too low
        wlo=whi; FElo=FEhi;
        whi=whi*1.3; FEhi=AV2001_FreeEntryResidual(whi,Params,Grids,vfoptions);
    end
    bracketiter=bracketiter+1;
end
if FElo*FEhi>0
    error('AV2001_StationaryEqm: could not bracket wstar for the free entry condition (i=%.6f)',i)
end
fzerooptions=optimset('TolX',1e-8*Params.wstar_guess,'Display','off');
wstar=fzero(@(ww) AV2001_FreeEntryResidual(ww,Params,Grids,vfoptions),[wlo,whi],fzerooptions);
[FEresidual,~,Policy_firm]=AV2001_FreeEntryResidual(wstar,Params,Grids,vfoptions);
Params.wstar=wstar;

% Employment decision n(e,s) as e-grid indices (CPU copy for our own indexing; the toolkit keeps its own) and as values
npol_idx=gather(shiftdim(Policy_firm(1,:,:),1)); % [n_e,n_s]
n_pol=e_grid(npol_idx); % [n_e,n_s]
if verbose>=2
    fprintf('    wstar=%.6f (free entry residual %.2e), max n index %d of %d \n',wstar,FEresidual,max(npol_idx,[],'all'),n_e)
end

%% Step 3: stationary distribution of establishments with entry mass nu=1
% Establishments at s=0 are in the distribution for one period (they fire
% everyone and pay the firing cost) and then disappear: survival prob is zero at s=0.
simoptions_firm=simoptions;
simoptions_firm.agententryandexit=1;
simoptions_firm.endogenousexit=0;
EntryExitParamNames.DistOfNewAgents={'firmentrydist'};
EntryExitParamNames.MassOfNewAgents={'nu'};
EntryExitParamNames.CondlProbOfSurvival={'firmsurvival'};
Params.firmentrydist=zeros(n_e,n_s);
Params.firmentrydist(1,2)=Params.psi(1);
Params.firmentrydist(1,3)=Params.psi(2);
Params.nu=1;
Params.firmsurvival=ones(n_e,n_s);
Params.firmsurvival(:,1)=0;
StationaryDist_firm=StationaryDist_InfHorz(Policy_firm,0,n_e,n_s,Q,simoptions_firm,Params,EntryExitParamNames);

FnsToEvaluate_firm.employment=@(n,e,s) n;
FnsToEvaluate_firm.capital=@(n,e,s,theta,gamma,r) (theta*s*n^gamma/r)^(1/(1-theta));
FnsToEvaluate_firm.output=@(n,e,s,theta,gamma,r) s*((theta*s*n^gamma/r)^(1/(1-theta)))^theta*n^gamma;
FnsToEvaluate_firm.hires=@(n,e,s,rho) max(n-(1-rho)*e,0);
FnsToEvaluate_firm.fires=@(n,e,s,rho) max((1-rho)*e-n,0);
AggVars_firm=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist_firm,Policy_firm,FnsToEvaluate_firm,Params,[],0,n_e,n_s,[],e_grid,Grids.s_grid,simoptions_firm,EntryExitParamNames);

firmdist=gather(StationaryDist_firm.mass*StationaryDist_firm.pdf); % [n_e,n_s], entry mass nu=1

%% Step 4: reachable establishment states, Lambda (probability of keeping the job) and G (distribution of new hires)
% Reachable set: entrants are at (e=0, s in {1,sbar}); from (e,s) the establishment moves to (n(e,s), s') for s' in {1,sbar}.
% (s=0 is death: everyone is fired, so it is never a household employed state.)
reach=false(n_e,n_s);
reach(1,2)=true; reach(1,3)=true;
frontier=[1,2;1,3];
while ~isempty(frontier)
    newfrontier=zeros(0,2);
    for kk=1:size(frontier,1)
        eprime_idx=npol_idx(frontier(kk,1),frontier(kk,2));
        for sprime_idx=2:3
            if ~reach(eprime_idx,sprime_idx)
                reach(eprime_idx,sprime_idx)=true;
                newfrontier(end+1,:)=[eprime_idx,sprime_idx]; %#ok<AGROW>
            end
        end
    end
    frontier=newfrontier;
end
[e_R,s_R]=find(reach); % the reachable (e_idx,s_idx) pairs, all s=1 first then all s=sbar
N_R=length(e_R);
hhidx=zeros(n_e,n_s); % (e_idx,s_idx) -> household employed state index
hhidx(sub2ind([n_e,n_s],e_R,s_R))=1:N_R;
offreachmass=sum(firmdist(~reach & repmat([false,true,true],n_e,1)),'all')/sum(firmdist(:,2:3),'all');
if offreachmass>1e-8
    warning('AV2001_StationaryEqm: establishment distribution has mass %.2e (share) on states not reachable from entry; the household state space misses them',offreachmass)
end

% Lambda, eq (3.4), as a function of the NEXT period establishment state (e',s'). If e'=n(e,s)=0 nobody can be kept.
Lambda=zeros(n_e,n_s);
Lambda(2:end,:)=min(1,n_pol(2:end,:)./((1-Params.rho)*e_grid(2:end)));
% G, eq (3.9), over the reachable states (hiring is zero everywhere else)
hiresmat=max(n_pol-(1-Params.rho)*e_grid,0).*firmdist; % [n_e,n_s]
G_R=hiresmat(sub2ind([n_e,n_s],e_R,s_R))/sum(hiresmat,'all'); % [N_R,1]
if abs(sum(G_R)-1)>1e-8
    error('AV2001_StationaryEqm: hiring distribution G over reachable states sums to %.6f, not one',sum(G_R))
end

%% Step 5: households at w=1
% Household state index: 1:N_R employed at reachable (e,s), then just-fired (severance in hand), unemployed, retired, dead.
% With unemployment insurance (Params.uirr>0, the calibrated U.S. economy of doPart(7)) the unemployed are also split by
% eligibility: just-fired ineligible, unemployed ineligible, just-fired eligible, unemployed eligible, retired, dead.
% Eligible agents receive uirr*w each period; a fired worker is eligible with probability uikappa, and an eligible
% unemployed agent who stays unemployed keeps eligibility with probability uiphi (once lost it is not regained).
UI=(Params.uirr>0);
if UI
    N_hh=N_R+6;
    z_fired=N_R+1; z_unemp=N_R+2; z_firedE=N_R+3; z_unempE=N_R+4; z_ret=N_R+5; z_dead=N_R+6;
else
    N_hh=N_R+4;
    z_fired=N_R+1; z_unemp=N_R+2; z_ret=N_R+3; z_dead=N_R+4;
end
semiz_grid=zeros(N_hh,1); % encoding 10*status+income, see AV2001_HouseholdReturnFn
semiz_grid(1:N_R)=1; % employed: wage
semiz_grid(z_fired)=10+Params.lambda; % unemployed, receives severance lambda*w this period
semiz_grid(z_unemp)=10;
if UI
    semiz_grid(z_firedE)=10+Params.lambda+Params.uirr; % severance and the first period of benefits
    semiz_grid(z_unempE)=10+Params.uirr;
end
semiz_grid(z_ret)=20;
semiz_grid(z_dead)=30;

% Rows that do not depend on search intensity
pibase=zeros(N_hh,N_hh);
for k=1:N_R
    e_idx=e_R(k); s_idx=s_R(k);
    eprime_idx=npol_idx(e_idx,s_idx); % next period the establishment has e'=n(e,s)
    pibase(k,z_ret)=Params.rho;
    pfired=(1-Params.rho)*Q(s_idx,1); % establishment dies: everyone fired
    for sprime_idx=2:3
        kprime=hhidx(eprime_idx,sprime_idx); % exists, the reachable set is closed under the policy
        keepprob=Lambda(eprime_idx,sprime_idx);
        pibase(k,kprime)=pibase(k,kprime)+(1-Params.rho)*Q(s_idx,sprime_idx)*keepprob;
        pfired=pfired+(1-Params.rho)*Q(s_idx,sprime_idx)*(1-keepprob);
    end
    if UI
        pibase(k,z_firedE)=Params.uikappa*pfired;
        pibase(k,z_fired)=(1-Params.uikappa)*pfired;
    else
        pibase(k,z_fired)=pfired;
    end
end
pibase_vfi=pibase; % value function: dead is absorbing with zero value, so retirees have no bequest motive
pibase_vfi(z_ret,z_ret)=Params.sigma;
pibase_vfi(z_ret,z_dead)=1-Params.sigma;
pibase_vfi(z_dead,z_dead)=1;
pibase_dist=pibase; % distribution: the offspring is born unemployed with the inherited assets
pibase_dist(z_ret,z_ret)=Params.sigma;
pibase_dist(z_ret,z_unemp)=1-Params.sigma;
pibase_dist(z_dead,z_unemp)=1; % no one ever dies in the distribution chain, but the toolkit's default initial distribution puts mass on every state, so the dead state must drain (if it were absorbing that mass would sit there forever, 1/N_hh of all agents)
% Unemployed rows (both with and without severance in hand) depend on eta through the job-finding probability eta^xi
pi_vfi=repmat(pibase_vfi,[1,1,n_eta]);
pi_dist=repmat(pibase_dist,[1,1,n_eta]);
for dd=1:n_eta
    jf=eta_grid(dd)^Params.xi;
    row=zeros(1,N_hh);
    row(z_ret)=Params.rho;
    row(1:N_R)=(1-Params.rho)*jf*G_R';
    row(z_unemp)=(1-Params.rho)*(1-jf);
    pi_vfi(z_fired,:,dd)=row; pi_vfi(z_unemp,:,dd)=row;
    pi_dist(z_fired,:,dd)=row; pi_dist(z_unemp,:,dd)=row;
    if UI % eligible: stay eligible with probability uiphi if still unemployed
        rowE=row;
        rowE(z_unemp)=(1-Params.rho)*(1-jf)*(1-Params.uiphi);
        rowE(z_unempE)=(1-Params.rho)*(1-jf)*Params.uiphi;
        pi_vfi(z_firedE,:,dd)=rowE; pi_vfi(z_unempE,:,dd)=rowE;
        pi_dist(z_firedE,:,dd)=rowE; pi_dist(z_unempE,:,dd)=rowE;
    end
end
if max(abs(sum(pi_vfi,2)-1),[],'all')>1e-10 || max(abs(sum(pi_dist,2)-1),[],'all')>1e-10
    error('AV2001_StationaryEqm: household transition matrix rows do not sum to one (max error %.2e)',max(abs(sum(pi_vfi,2)-1),[],'all'))
end
pi_vfi=pi_vfi./sum(pi_vfi,2);
pi_dist=pi_dist./sum(pi_dist,2);

Params.w=1;
if Params.endogsearch==1
    if verbose>=1
        N_aprimefine=(1+2*Grids.maxaprimediff)+2*Grids.maxaprimediff*Grids.ngridinterp;
        fprintf('    household: N_R=%d reachable establishment states, N_hh=%d, coarse return matrix %.2f GB, each fine array %.2f GB \n',N_R,N_hh,n_eta*n_a*n_a*N_hh*8/1e9,N_aprimefine*n_a*N_hh*n_eta*8/1e9*(Grids.ngridinterp>0))
    end
    vfoptions_hh=vfoptions;
    % Howards by indexing rather than with a sparse matrix once the household has many states: the toolkit builds the
    % sparse matrix from triplets with the transition matrix expanded densely, N_a*N_hh*N_hh entries (and twice that with
    % the interpolation layer), although each household state moves to at most five others. With 860 states (the
    % attrition-exact employment grid of doPart(8d)) that is ~25-30 GB and ran out of memory (2026-10-28); indexing
    % peaks at ~5 GB there.
    vfoptions_hh.howardssparse=(N_hh<=Grids.howardssparse_maxNhh);
    vfoptions_hh.n_semiz=N_hh;
    vfoptions_hh.semiz_grid=semiz_grid;
    vfoptions_hh.pi_semiz=pi_vfi;
    simoptions_hh=simoptions;
    simoptions_hh.agententryandexit=0;
    simoptions_hh.n_semiz=N_hh;
    simoptions_hh.semiz_grid=semiz_grid;
    simoptions_hh.pi_semiz=pi_dist;
    simoptions_hh.d_grid=eta_grid;
    if Grids.ngridinterp>0 % grid interpolation layer on assets (the coarse grid is far too coarse for the per-period saving flows of a half-quarter model); maxaprimediff must be set: the toolkit default is a fifth of the grid, which with a decision variable is a 60+ GB fine array
        vfoptions_hh.gridinterplayer=1; vfoptions_hh.ngridinterp=Grids.ngridinterp; vfoptions_hh.maxaprimediff=Grids.maxaprimediff; vfoptions_hh.verbose_advice=0;
        simoptions_hh.gridinterplayer=1; simoptions_hh.ngridinterp=Grids.ngridinterp;
    end
    HouseholdReturnFn=@(eta,aprime,a,z,i,w,alpha,tau) AV2001_HouseholdReturnFn_endog(eta,aprime,a,z,i,w,alpha,tau);
    [V_hh,Policy_hh]=ValueFnIter_InfHorz(n_eta,n_a,0,eta_grid,a_grid,[],[],HouseholdReturnFn,Params,{'beta'},[],vfoptions_hh);
    StationaryDist_hh=StationaryDist_InfHorz(Policy_hh,n_eta,n_a,0,[],simoptions_hh);

    FnsToEvaluate_hh.consumption=@(eta,aprime,a,z,i,w) (1+i)*a+w*(z-10*floor(z/10))-aprime;
    FnsToEvaluate_hh.assets=@(eta,aprime,a,z) a;
    FnsToEvaluate_hh.employed=@(eta,aprime,a,z) (floor(z/10)==0);
    FnsToEvaluate_hh.unemployed=@(eta,aprime,a,z) (floor(z/10)==1);
    FnsToEvaluate_hh.retired=@(eta,aprime,a,z) (floor(z/10)==2);
    FnsToEvaluate_hh.searchtime=@(eta,aprime,a,z) eta*(floor(z/10)==1);
    FnsToEvaluate_hh.jobfinders=@(eta,aprime,a,z,xi) (eta^xi)*(floor(z/10)==1);
    AggVars_hh=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist_hh,Policy_hh,FnsToEvaluate_hh,Params,[],n_eta,n_a,0,eta_grid,a_grid,[],simoptions_hh);
    FnsToEvaluate_hh2.consumption=FnsToEvaluate_hh.consumption;
    AllStats_hh=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist_hh,Policy_hh,FnsToEvaluate_hh2,Params,[],n_eta,n_a,0,eta_grid,a_grid,[],simoptions_hh);
    W_hh_noalpha=NaN; % with endogenous search, alpha=0 changes the policy; the paper only reports this row for the exogenous economies
else
    pi_z_vfi=pi_vfi(:,:,1);
    pi_z_dist=pi_dist(:,:,1);
    vfoptions_hh=vfoptions;
    vfoptions_hh.howardssparse=(N_hh<=Grids.howardssparse_maxNhh); % see the endogenous branch
    simoptions_hh=simoptions;
    simoptions_hh.agententryandexit=0;
    if Grids.ngridinterp>0 % grid interpolation layer on assets, see the endogenous branch
        vfoptions_hh.gridinterplayer=1; vfoptions_hh.ngridinterp=Grids.ngridinterp; vfoptions_hh.maxaprimediff=Grids.maxaprimediff; vfoptions_hh.verbose_advice=0;
        simoptions_hh.gridinterplayer=1; simoptions_hh.ngridinterp=Grids.ngridinterp;
    end
    HouseholdReturnFn=@(aprime,a,z,i,w,alpha,tau,eta) AV2001_HouseholdReturnFn(aprime,a,z,i,w,alpha,tau,eta);
    [V_hh,Policy_hh]=ValueFnIter_InfHorz(0,n_a,N_hh,[],a_grid,semiz_grid,pi_z_vfi,HouseholdReturnFn,Params,{'beta'},[],vfoptions_hh);
    StationaryDist_hh=StationaryDist_InfHorz(Policy_hh,0,n_a,N_hh,pi_z_dist,simoptions_hh);

    FnsToEvaluate_hh.consumption=@(aprime,a,z,i,w) (1+i)*a+w*(z-10*floor(z/10))-aprime;
    FnsToEvaluate_hh.assets=@(aprime,a,z) a;
    FnsToEvaluate_hh.employed=@(aprime,a,z) (floor(z/10)==0);
    FnsToEvaluate_hh.unemployed=@(aprime,a,z) (floor(z/10)==1);
    FnsToEvaluate_hh.retired=@(aprime,a,z) (floor(z/10)==2);
    FnsToEvaluate_hh.searchtime=@(aprime,a,z,eta) eta*(floor(z/10)==1);
    FnsToEvaluate_hh.jobfinders=@(aprime,a,z,eta,xi) (eta^xi)*(floor(z/10)==1);
    AggVars_hh=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist_hh,Policy_hh,FnsToEvaluate_hh,Params,[],0,n_a,N_hh,[],a_grid,semiz_grid,simoptions_hh);
    FnsToEvaluate_hh2.consumption=FnsToEvaluate_hh.consumption;
    AllStats_hh=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist_hh,Policy_hh,FnsToEvaluate_hh2,Params,[],0,n_a,N_hh,[],a_grid,semiz_grid,simoptions_hh);

    % Welfare with no search disutility: with eta exogenous, alpha only adds a status-dependent constant to V, so the
    % policy is unchanged and re-solving the value function with alpha=0 is exact.
    Params_noalpha=Params; Params_noalpha.alpha=0;
    V_hh_noalpha=ValueFnIter_InfHorz(0,n_a,N_hh,[],a_grid,semiz_grid,pi_z_vfi,HouseholdReturnFn,Params_noalpha,{'beta'},[],vfoptions_hh);
    W_hh_noalpha=sum(gather(V_hh_noalpha).*gather(StationaryDist_hh),'all');
end

% Average lifetime utility across agents (w=1 units)
W_hh=sum(gather(V_hh).*gather(StationaryDist_hh),'all');

N_hh_emp=AggVars_hh.employed.Mean;
U_hh=AggVars_hh.unemployed.Mean;
R_hh=AggVars_hh.retired.Mean;
jobfinders_mass=AggVars_hh.jobfinders.Mean; % unemployed who find the labour market (before retirement is applied)
jobfindingrate=jobfinders_mass/U_hh; % average eta^xi across the unemployed
nonsearchtime=1-AggVars_hh.searchtime.Mean;

%% Step 6: scale of the establishment distribution nu, from labour market clearing, eq (A.1)
% Search decisions are homogeneous of degree zero in w, so the household side at w=1 is already correct.
jobfinders=(1-Params.rho)*jobfinders_mass;
nu=jobfinders/sum(hiresmat,'all');

%% Step 7: aggregates of establishments at the correct scale
N_firm=nu*AggVars_firm.employment.Aggregate;
K=nu*AggVars_firm.capital.Aggregate;
Y=nu*AggVars_firm.output.Aggregate;
Hires=nu*AggVars_firm.hires.Aggregate;
Fires=nu*AggVars_firm.fires.Aggregate;
massfirms=nu*StationaryDist_firm.mass;

%% Step 8: after payroll tax wage w from the government budget constraint, eq (3.13)
% (wstar-w)*N + (epsilon*wstar-lambda*w)*Fires - uirr*w*U_eligible = 0
% (with unemployment insurance the payroll tax also pays the benefits the firing taxes do not)
if UI
    statemass=sum(gather(StationaryDist_hh),1); % [1,N_hh], mass in each household state
    U_elig=statemass(z_firedE)+statemass(z_unempE); % agents receiving benefits this period
    UIduration=U_elig/statemass(z_firedE); % periods of benefits per spell of eligibility (stock over inflow)
else
    U_elig=0; UIduration=0;
end
w=wstar*(N_hh_emp+Params.epsilon*Fires)/(N_hh_emp+Params.lambda*Fires+Params.uirr*U_elig);

%% Step 9: consumption is homogeneous of degree one in w
C=w*AggVars_hh.consumption.Mean;
S=w*AggVars_hh.assets.Mean;
sdC=w*AllStats_hh.consumption.StdDeviation;

%% Step 10: goods market clearing, eq (3.14)
residual=(C+Params.delta*K+Params.m*nu-Y)/Y;

%% Welfare: log utility, so V at wage w is V(w=1)+log(w)*L(status), where L is the
% expected discounted lifetime (actives and retirees differ because of retirement and death).
L_R=1/(1-Params.beta*Params.sigma);
L_A=(1+Params.beta*Params.rho*L_R)/(1-Params.beta*(1-Params.rho));
Lbar=(N_hh_emp+U_hh)*L_A+R_hh*L_R;
W=W_hh+log(w)*Lbar;
W_noalpha=W_hh_noalpha+log(w)*Lbar;

%% Collect
Out.i=i; Out.annualinterestrate=(1+i)^8-1;
Out.wstar=wstar; Out.w=w; Out.nu=nu; Out.massfirms=massfirms;
Out.avgestablishmentsize=N_firm/massfirms;
Out.jobfindingrate=jobfindingrate;
Out.layoffrate=Fires/N_hh_emp;
Out.unemploymentrate=U_hh/(N_hh_emp+U_hh);
Out.nonsearchtime=nonsearchtime;
Out.N_hh=N_hh_emp; Out.N_firm=N_firm; Out.U=U_hh; Out.R=R_hh;
Out.U_elig=U_elig; Out.UIduration=UIduration; Out.UIweeks=6.5*UIduration; % a model period is half a quarter, 6.5 weeks
Out.firingtaxshareofUI=(Params.uirr>0)*Params.epsilon*wstar*Fires/max(Params.uirr*w*U_elig,eps); % share of UI benefits paid by firing taxes
Out.Hires=Hires; Out.Fires=Fires;
Out.Y=Y; Out.K=K; Out.C=C; Out.S=S; Out.sdC=sdC;
% Measured output excludes the investment in new establishments (Section 4: "the entrepreneurial investment in new
% establishments goes unmeasured"), as does the calibration target K/Y=1.7.
Out.Ymeasured=Y-Params.m*nu;
Out.KYannual=K/(8*Out.Ymeasured);
% Value of the assets owned by banks, eq (3.19): capital plus the value of the establishments, the present value of
% their profits net of entry. In equilibrium A=S (goods market clearing implies it).
profits=Y-wstar*N_firm-Params.r*K-Params.epsilon*wstar*Fires;
Out.A=K+(profits-nu*Params.m)/i;
Out.W=W; Out.W_noalpha=W_noalpha; Out.Lbar=Lbar;
Out.residual=residual;
Out.FEresidual=FEresidual;
Out.N_R=N_R; Out.N_hhstates=N_hh;
Out.masstopassetgrid=sum(gather(StationaryDist_hh(end,:)),'all');
Out.masstopegrid=sum(firmdist(end,:),'all')/sum(firmdist,'all');
Out.n_pol=n_pol; Out.firmdist=nu*firmdist; Out.reach=reach;

if verbose>=1
    fprintf('  i=%.6f (annual %.3f%%): wstar=%.5f w=%.5f nu=%.4f, resid=%.3e, jf=%.4f, u=%.4f, layoff=%.4f, K/Y(annual)=%.3f, avg size=%.1f, N_hh=%.4f N_firm=%.4f, N_R=%d, mass top a-grid=%.2e, mass top e-grid=%.2e \n', ...
        i,100*Out.annualinterestrate,wstar,w,nu,residual,jobfindingrate,Out.unemploymentrate,Out.layoffrate,Out.KYannual,Out.avgestablishmentsize,N_hh_emp,N_firm,N_R,Out.masstopassetgrid,Out.masstopegrid)
end

end
