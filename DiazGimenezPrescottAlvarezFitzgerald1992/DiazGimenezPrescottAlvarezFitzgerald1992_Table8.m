function Table8=DiazGimenezPrescottAlvarezFitzgerald1992_Table8;
% For Table 8 the model period is now one-year (rather than one-eigth of a year)
Table8=zeros(3,7);

% Agents no longer die, they live forever

n_s=2;
n_z=1;

%Real wage, w(s,z)
Params.w1=1;
Params.w2=0.4;

% Parameter values from pg 555
%Calibrate household parameters
%Preferences
Params.alpha=0.3330; %Private consumption share
Params.beta=0.96;
Params.psi=4.0000; %Risk aversion
Params.tau=2.2200; %Productive time
Params.delta_g=0.039111; %Public consumption factor

Params.rho=1;     %reserve requirement(z)
Params.theta=0.2; %Tax rate on labour & interest income(z)

% Baseline zero inflation
Params.epsilon_inflation=1.000;    %Inflation rate process

% The real return on deposits is equal to minus the inflation rate:
Params.eta_d=0; %Deposits (to match Imrohoroglu & Prescott 1991)

Params.iota=0.00625;    %Nominal interest rate on T-bill(z) 

%w the Real wage, w(s,z) is declared above

%calculate a few of the other parameters directly
Params.i_d=(1-Params.rho)*Params.iota-Params.eta_d;%=0 in practice anyway %interest rate on deposits
Params.e=Params.epsilon_inflation-1; %pricing/inflation process on reserves



% Following parameters are irrelevant as no longer allowed loans nor
% capital; nor is there retirement nor death.
Params.alpha_k=0.1080; %Capital service share
Params.beta=0.9994; %Time discount factor
Params.delta_r=0.2100; %Retirees' constant
%Technology
Params.mew=0.00625; %Maintainence cost                     %%%%Note: Value in table in paper is incorrect. Paper reports annual, not eighth-of-year, value.
Params.gamma=0.02625; %Rental service coefficient       %%%% Pg., 550: "we set gamma to twice the sum of the real borrowing rate (i_l-e) and the maintainence cost (mew)"
Params.phi=0.9000; %1-phi = 0.1000; Disinvestment cost
Params.phi_1=0.9000; %Prob of newborn being: Type 1
Params.phi_2=0.1000; %                       Type 2
%Per unit banking costs
Params.eta_l=0.0056250; %Loans
%Government policy declared above
%Calibrated bank and government parameters
Params.omega=0.02; % Welfare tranfers to Indigent retirees (are zero to everyone else)
Params.i_l=Params.iota+Params.eta_l; %interest rate on loans
Params.w3=0; Params.w4=0; Params.Experiment=0;

%% Create grids for s,z,a (A & K), and transition matrices for s & z
n_A=800; %number taken from tech append II, pg 26
% PROHIBIT Borrowing from banks by setting lower bound on A_grid to be zero.
% Other than restricting to no negative assests (so as no capital there is
% no lending) I leave the maximum value of A unchanged, and number of
% points unchanged. This means that grid will be different.
A_grid=linspace(0,(2.7/245)*(n_A-1)-2.7,n_A)'; %so we n_A points
[~,n_A_zero]=min(abs(A_grid-0)); % find index that corresponds to zero assets (the minus zero is redundant here, but used to make exposition clearer)
K_grid=0; %household capital stocks (pg 550) (house + consumer durables + small business capital)
n_K=1;

N_grid=[0;1]; %to use for calculating the labour from it's index


%Calibrated household idiosyncratic transtion probablities
pi_s=[0.9565, 1-0.9565; 0.5000, 0.5000]; % pg 555
s_grid=(1:1:2)';  % 2 individual states
% transition matrix for s, the individual specific states row index is this state s, column index is next state s'
%Economy-wide shock transition probabilities
pi_z=1;
z_grid=1;
w_sz=[Params.w1;Params.w2]; % used for welfare evaluations
sigma_sz=[1;1]; % used for welfare evaluations

%Note that from the point of view of the value function, there is no difference between s & z, so we combine them.
sz_grid=[s_grid;z_grid];
pi_sz=kron(pi_s,pi_z);
n_s=length(s_grid);
n_z=length(z_grid);
n_sz=[n_s,n_z];

a_grid=[A_grid;K_grid];
n_a=[length(A_grid),length(K_grid)];

d_grid=N_grid;
n_d=length(N_grid);

DiscountFactorParamNames={'beta'};

ReturnFn=@(N_val,Aprime_val,Kprime_val, A_val,K_val, s_val,z_val,phi,e,omega,w1,w2,w3,w4,theta,i_d,i_l,mew,alpha,alpha_k,gamma,tau,psi,delta_r,Experiment) DiazGimenezPrescottAlvarezFitzgerald1992_ReturnFn(N_val,Aprime_val,Kprime_val, A_val,K_val, s_val,z_val,phi,e,omega,w1,w2,w3,w4,theta,i_d,i_l,mew,alpha,alpha_k,gamma,tau,psi,delta_r,Experiment)
ReturnFnParamNames={'phi','e','omega','w1','w2','w3','w4','theta','i_d','i_l','mew','alpha','alpha_k','gamma','tau','psi','delta_r','Experiment'};


%%

% Set inflation rates to get real after-tax returns of 0,-2,-4 & -6 percent.
% Since pre-tax return on deposits is zero (as 100% reserves, rho=1, and no
% intermediation cost, eta_d=0) the after-tax nominal return will just be
% equal to the zero pre-tax nominal return, and so real after-tax return is
% just the negative of the inflation rate.
inflation_vec=[0,0.02,0.04,0.06];


% Generate objects needed to calculate rows of Table 8 (the reference case and the three rows)
for ii=1:4 
    Params.epsilon_inflation=1+inflation_vec(ii);    %Inflation rate process
    
    Params.i_d=(1-Params.rho)*Params.iota-Params.eta_d;%=0 in practice anyway %interest rate on deposits
    Params.e=Params.epsilon_inflation-1; %pricing/inflation process on reserves
    
    V0=ones([n_a,n_sz]);
    vfoptions.policy_forceintegertype=1;
    [V, Policy]=ValueFnIter_Case1(V0, n_d,n_a,n_sz,d_grid,a_grid,sz_grid, pi_sz, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames,vfoptions); %
    
    % Since there is no longer the death and rebirth can just use standard VFI
    % Toolkit codes for getting agent distribution.
    
    StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_sz,pi_sz);
    
    [ValueFnofPublicandPrivateConsumption,v1,v2,W]=DiazGimenezPrescottAlvarezFitzgerald1992_Welfare(Params,V,n_A,n_K,n_s,n_z,n_d,n_a,n_sz,sigma_sz,w_sz,pi_sz,A_grid,K_grid,d_grid, a_grid,sz_grid,StationaryDist, Policy);

    Outputs(ii).V=gather(V);
    Outputs(ii).Policy=gather(Policy);
    Outputs(ii).StationaryDist=gather(StationaryDist);
    Outputs(ii).ValueFnofPublicandPrivateConsumption=gather(ValueFnofPublicandPrivateConsumption);
    Outputs(ii).v1=gather(v1);
    Outputs(ii).v2=gather(v2);
    Outputs(ii).W=gather(W);
    
    % Compute the expected utility across agents (mean of value function at
    % stationary distribution). This is the basis of the
    % consumption-equivalence measures used by Imrohoroglu & Prescott
    % (1991) that form the left-half of Table 8.
    Outputs(ii).EV=sum(sum(sum(Outputs(ii).V.*Outputs(ii).StationaryDist)));
    Outputs(ii).EValueFnofPublicandPrivateConsumption=sum(sum(sum(Outputs(ii).ValueFnofPublicandPrivateConsumption.*Outputs(ii).StationaryDist)));
    Outputs(ii).Ev1=sum(sum(sum(Outputs(ii).v1.*Outputs(ii).StationaryDist)));
    Outputs(ii).Ev2=sum(sum(sum((ones(n_A,n_K).*shiftdim(Outputs(ii).v2,-2)).*Outputs(ii).StationaryDist)));
    Outputs(ii).EW=sum(sum(sum(Outputs(ii).W.*Outputs(ii).StationaryDist)));
end

% Now calculate the actual entries
for ii=2:4
    % First, calculate the measure used by the paper of evaluating at each
    % point of the actual stationary distribution, and then taking mean
    % across this.
    % [Paper incorrectly describes this as taking account of the
    % 'transitional dynamics'. It is not accounting for transition path.
    % It accounts for differences in agents across the stationary distribution
    % of agents.]
    Lambda=(Outputs(1).ValueFnofPublicandPrivateConsumption ./ Outputs(ii).ValueFnofPublicandPrivateConsumption).^(1/(Params.alpha*(1-Params.psi)));
    % n_z=1, so can ignore this aspect of formula for M (pg 554) (y in
    % formula on pg 554 is the stationary distribution for the pi0 policy)
    M=sum(sum(sum(Outputs(ii).W.*(Lambda-1).*Outputs(ii).StationaryDist)))/sum(sum(sum(Outputs(ii).W.*Outputs(ii).StationaryDist)));
    
    M1factor=(Outputs(1).v1-Outputs(ii).v1)./(Outputs(1).ValueFnofPublicandPrivateConsumption-Outputs(ii).ValueFnofPublicandPrivateConsumption);
    M2factor=(ones(n_A,n_K).*shiftdim(Outputs(1).v2-Outputs(ii).v2,-2))./(Outputs(1).ValueFnofPublicandPrivateConsumption-Outputs(ii).ValueFnofPublicandPrivateConsumption);
    
    M1=sum(sum(sum(Outputs(ii).W.*M1factor.*(Lambda-1).*Outputs(ii).StationaryDist)))/sum(sum(sum(Outputs(ii).W.*Outputs(ii).StationaryDist)));
    M2=sum(sum(sum(Outputs(ii).W.*M2factor.*(Lambda-1).*Outputs(ii).StationaryDist)))/sum(sum(sum(Outputs(ii).W.*Outputs(ii).StationaryDist)));

    Table8(ii-1,5)=M;
    Table8(ii-1,6)=M1;
    Table8(ii-1,7)=M2;
    
    
    % Second, calculate in the way used by Prescott & Imrohorglu (1991) of
    % taking means, and then performing the calculation to find the
    % compensating change.
    % [Paper incorrectly describes this as the 'steady-state' when it is in
    % fact valued at mean of stationary distribution.]
    LambdaE=(Outputs(1).EValueFnofPublicandPrivateConsumption./Outputs(ii).EValueFnofPublicandPrivateConsumption).^(1/(Params.alpha*(1-Params.psi)));
    S=(Outputs(ii).EW*(LambdaE-1))/Outputs(ii).EW; %Obviously this just gives LambdaE-1, I keep this form to highlight how S relates to M.
    
    S1factor=(Outputs(1).Ev1-Outputs(ii).Ev1)./(Outputs(1).EValueFnofPublicandPrivateConsumption-Outputs(ii).EValueFnofPublicandPrivateConsumption);
    S2factor=(Outputs(1).Ev2-Outputs(ii).Ev2)./(Outputs(1).EValueFnofPublicandPrivateConsumption-Outputs(ii).EValueFnofPublicandPrivateConsumption);
    
    S1=Outputs(ii).EW.*S1factor.*(LambdaE-1)/Outputs(ii).EW;
    S2=Outputs(ii).EW.*S2factor.*(LambdaE-1)/Outputs(ii).EW;
    
    Table8(ii-1,2)=S;
    Table8(ii-1,3)=S1;
    Table8(ii-1,4)=S2;
end


Table8(:,1)=-inflation_vec(2:4)';


Table8=100*Table8;

% Some toying around I was doing related to results of Prescott & Imrohoroglu (1991) and their results
% u0=-0.2847
% u1=-0.2863
% psi=4; alpha=0.333;
% lamb=(u1/u0)^(1/(alpha*(1-psi)))
% M_Expected=lamb-1








