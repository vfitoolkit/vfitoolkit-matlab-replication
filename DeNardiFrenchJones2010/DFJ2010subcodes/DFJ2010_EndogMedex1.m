% DFJ2010_EndogMedex1: the model with ENDOGENOUS medical expenses (DFJ2010 Sec. IX, eqns 21-24; online appendix A.2, B, C;
% replication materials singlsav51_wave8.gau, gmmcrit51.g, getdat50_wave8.g, wealth45_Singles.c; see
% EndogenousMedex_proposal.md), part 1: parameters at DFJ2010's published endogenous-model estimates (Table A1 column 5
% and the medical-needs shifter in singlsav51_wave8.gau), grids and shock transitions, return function, value function.
% Called from DeNardiFrenchJones2010.m as doPart(6); expects the first-stage parameters (DFJ2010_FirstStageParams) in
% the workspace. OVERWRITES the exogenous model's Params.rho_m/sigma_mepsilon/sigma_mxi/normalizepsi, z_grid, pi_z_J,
% ReturnFn, V and Policy (this part runs after parts 2-5, which are unaffected).
%
% The model keeps the state space of the exogenous model (assets; health h and the persistent shock zeta as markov
% states; the iid shock xi) and the single choice a'. Medical goods m are chosen too, but the intratemporal first-order
% condition m=(mu/q)^(1/omega)*c^(nu/omega) makes the split a scalar root, found by bisection inside the return function
% (as DFJ2010's C code does), so no decision variable is needed.

% Preferences: DFJ2010 online appendix Table A1 column 5 (singlsav51_wave8.gau lines 217, 324-326)
Params.nu=2.14608;      % 1/IES for consumption (paper: 2.15, se .055)
Params.omega=3.19104;   % 1/IES for medical goods (paper: 3.19, se .054)
Params.beta=0.990474;   % discount factor (paper: .99)
Params.eMin=exp(5.30863); % =202: the consumption level that delivers the utility floor when mu=1 (Table A1 footnote); the floor
                          % utility is uMin=u(eMin,eMin,mu=1) (C: uMin=getUtility(eMin,eMin,1)); paper: 202 (se 23.4)
% No bequest motive (phi0=1) and no health preference shifter (delta=0) in the estimated endogenous model.
DiscountFactorParamNames={'beta'};

%% The medical-needs shifter mu(t,h,zeta,xi) and its shocks (singlsav51_wave8.gau lines 340-365: mxvec)
% ln mu = D(t,h,male) + S, S = zeta + xi. mxvec(1:2) are the logit-transformed autocorrelation of zeta and the share of
% the (unit) variance of S that is persistent (gmmcrit51.g punscale_m); mxvec(3:13) the coefficients of the mean of
% ln mu on [1, age, age^2/100, age^3/10000, bad health, bad*age, male, male*age, PI, PI*age, PI^2] (getdat50_wave8.g
% getmxtab, age in years); mxvec(14:24) the same for the VARIANCE of S (the C code multiplies zeta+xi by its square root).
% The PI coefficients are fixed at zero (the paper: 'we do not allow any of these shocks to depend on permanent income').
mxvec=[1.72357, -0.0478949, -11.1505, 0.26789, -0.262621, 0.221613, -0.446767, 0.00339203, 0, 0, 0, 0, 0, ...
       21.5593, -0.132871, 0.238177, 0, -0.740969, 0.0203848, 0, 0, 0, 0, 0];
logit=@(x) exp(x)./(1+exp(x)); % getdat50_wave8.g line 5
Params.rho_m=logit(mxvec(1));          % 0.8486 (exogenous model: 0.922)
fracar1=logit(mxvec(2));               % 0.4880: share of var(S) from zeta (exogenous model: 0.335)
Params.sigma_mepsilon=sqrt(fracar1*(1-Params.rho_m^2)); % innovation std of zeta (gmmcrit25.g line 1268: fracar1i=fracar1*(1-rho^2))
Params.sigma_mxi=sqrt(1-fracar1);      % std of xi
Params.normalizepsi=sqrt(Params.sigma_mepsilon^2/(1-Params.rho_m^2)+Params.sigma_mxi^2); % =1 by construction
fprintf('Endogenous medical-needs shocks: rho %.4f, persistent share %.4f, sigma_eps %.4f, sigma_xi %.4f, var(psi) %.4f\n',Params.rho_m,fracar1,Params.sigma_mepsilon,Params.sigma_mxi,Params.normalizepsi^2)
ageseq=(70:102); % actual age (getdat50: ageseq=seqa(bornage,1,_tr))
for ii=1:length(Names_i)
    male=strncmp(Names_i{ii},'male',4);
    Dmean_good=mxvec(3)+mxvec(4)*ageseq+mxvec(5)*ageseq.^2/100+mxvec(6)*ageseq.^3/10000+male*(mxvec(9)+mxvec(10)*ageseq);
    Dmean_bad=Dmean_good+mxvec(7)+mxvec(8)*ageseq;
    Var_good=mxvec(14)+mxvec(15)*ageseq+mxvec(16)*ageseq.^2/100+mxvec(17)*ageseq.^3/10000+male*(mxvec(20)+mxvec(21)*ageseq);
    Var_bad=Var_good+mxvec(18)+mxvec(19)*ageseq;
    Params.mu_coeff_healthgood.(Names_i{ii})=Dmean_good;
    Params.mu_coeff_healthbad.(Names_i{ii})=Dmean_bad;
    Params.s_coeff_healthgood.(Names_i{ii})=sqrt(abs(Var_good)); % getmxtab: stdmxs=sqrt(mxvarscl*abs(var)), mxvarscl=1
    Params.s_coeff_healthbad.(Names_i{ii})=sqrt(abs(Var_bad));
end
fprintf('Medical-needs shifter, female good health: mean of ln mu at 74/90/100 = %s, std of ln mu = %s (E[mu] at 74 = %.0f)\n',mat2str(round(Params.mu_coeff_healthgood.femaleq1([74,90,100]-69),2)),mat2str(round(Params.s_coeff_healthgood.femaleq1([74,90,100]-69),2)),exp(Params.mu_coeff_healthgood.femaleq1(5)+0.5*Params.s_coeff_healthgood.femaleq1(5)^2))

%% Coinsurance (out-of-pocket share) q(t,h,male): data/medpriceprof.out, columns age, b_const, b_badhealth, b_male;
% q=exp(b_const+b_badhealth*bad+b_male*male) (getdat50_wave8.g lines 192-211, appendix B: 'the coinsurance rate falls
% with age and for people in bad health')
% Contents of medpriceprof.out (ages 70..102), pasted so the run does not need the replication materials present:
qprof.age=70:1:102;
qprof.b_const=[-1.094554, -1.112593, -1.134112, -1.157964, -1.183139, -1.208691, -1.233812, -1.257761, -1.279934, -1.299792, -1.316934, -1.331026, -1.341869, -1.349332, -1.35342, -1.354205, -1.351897, -1.346769, -1.339234, -1.32977, -1.318993, -1.307585, -1.296363, -1.286214, -1.278158, -1.273285, -1.272819, -1.278052, -1.290412, -1.311396, -1.342633, -1.385824, -1.442802];
qprof.b_badhealth=[-.6711934, -.6817632, -.6923331, -.7029029, -.7134728, -.7240426, -.7346125, -.7451823, -.7557521, -.766322, -.7768918, -.7874617, -.7980315, -.8086014, -.8191712, -.8297411, -.8403109, -.8508807, -.8614506, -.8720204, -.8825902, -.8931601, -.9037299, -.9142998, -.9248696, -.9354395, -.9460093, -.9565791, -.967149, -.9777188, -.9882887, -.9988585, -1.009428];
qprof.b_male=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
for ii=1:length(Names_i)
    male=strncmp(Names_i{ii},'male',4);
    Params.q_healthgood.(Names_i{ii})=exp(qprof.b_const+male*qprof.b_male);
    Params.q_healthbad.(Names_i{ii})=exp(qprof.b_const+qprof.b_badhealth+male*qprof.b_male);
end
fprintf('Coinsurance rate q, female: good health %s, bad health %s at ages 74/90/100\n',mat2str(round(Params.q_healthgood.femaleq1([74,90,100]-69),3)),mat2str(round(Params.q_healthbad.femaleq1([74,90,100]-69),3)))

%% Grids and shock transitions (shared subcode; uses the endogenous rho_m, sigma_mepsilon, sigma_mxi set above)
DFJ2010_Grids

%% ReturnFn
ReturnFn=@(aprime,a,h,zeta,xi,r,nu,omega,eMin,earnings,mu_coeff_healthbad,mu_coeff_healthgood,s_coeff_healthbad,s_coeff_healthgood,normalizepsi,q_healthbad,q_healthgood,taxbracket1,taxbracket2,taxbracket3,taxbracket4,taxbracket5,taxbracket6,margtaxrate0,margtaxrate1,margtaxrate2,margtaxrate3,margtaxrate4,margtaxrate5,margtaxrate6)...
    DFJ2010EM_ReturnFn(aprime,a,h,zeta,xi,r,nu,omega,eMin,earnings,mu_coeff_healthbad,mu_coeff_healthgood,s_coeff_healthbad,s_coeff_healthgood,normalizepsi,q_healthbad,q_healthgood,taxbracket1,taxbracket2,taxbracket3,taxbracket4,taxbracket5,taxbracket6,margtaxrate0,margtaxrate1,margtaxrate2,margtaxrate3,margtaxrate4,margtaxrate5,margtaxrate6);

%% Solve the endogenous model at the published parameters
vfoptions.verbose=1;
vfoptions.verboseparams=0;
vfoptions.divideandconquer=1;
vfoptions.gridinterplayer=1;
vfoptions.ngridinterp=10;
simoptions.gridinterplayer=vfoptions.gridinterplayer;
simoptions.ngridinterp=vfoptions.ngridinterp;
tic;
[V,Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
vftime=toc;
fprintf('Value fn (endogenous medical expenses, 10 ptypes) took %.1f seconds (exogenous model: about 20 s)\n',vftime)

