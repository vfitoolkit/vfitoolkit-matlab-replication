function Output=GuerrieriLorenzoni2017_Fn(n_d,n_a,n_theta,n_r,Params, vfoptions, simoptions)
%This code replicates the results of Guerrieri & Lorenzoni (2017) - Credit Crises, Precautionary Savings, and the Liquidity Trap
disp('Running GuerrieriLorenzoni2017_Fn')
%% Set up

%Create markov process for the exogenous income (based on idea of employment and unemployment states, following Imrohoroglu, 1989).
[theta1_grid, pi_theta1]=TauchenMethod(0, Params.sigma_epsilon.^2,Params.rho, n_theta-1, Params.tauchenq);
z_grid=[0; exp(theta1_grid)];
pistar_theta1=ones(n_theta-1,1)/(n_theta-1);
for ii=1:10^4 % G&L2017, pg 1438 "when first employed, workers draw theta from its unconditional distribution"
    pistar_theta1=pi_theta1'*pistar_theta1; % There is a more efficient form to do this directly from a formula but I am feeling lazy. %FIX THIS LATER!!!
end
pi_z=[(1-Params.pi_ue), Params.pi_ue*pistar_theta1'; Params.pi_eu*ones(n_theta-1,1),(1-Params.pi_eu)*pi_theta1];
% Rows were did not sum to one due to rounding errors at order of 10^(-11), fix this
pi_z=pi_z./sum(pi_z,2);
pistar_z=ones(n_theta,1)/n_theta;
for ii=1:10^4 %  % There is a more efficient way to do this directly from a formula but I am feeling lazy. %FIX THIS LATER!!!
    pistar_z=pi_z'*pistar_z; % Formula could be used to find stationary dist of the employment unemployment process, then just combine with stationary dist of theta1, which is already calculated
end
% "The average level of theta is chosen so that yearly output in the initial steady state is normalized to 1"
z_grid=z_grid/sum(z_grid.*pistar_z);
% Double-check that this is 1
% sum(z_grid.*pistar_z)

% That the "normalized to 1" refers to E[theta] and not E[n*theta] is clear from setting
% v=0.1 to satisfy "For the unemployment benefit, we also follow Shimer
% (2005) and set it to 40% of average labor income." (pg 1438)

%% Grids
% Set grid for asset holdings
Params.alowerbar=-2; % This seems reasonable (No-one can go below -Params.phi in any case)
Params.aupperbar=12; % Not clear exactly what value is appropriate, have gone with this based on axes of Figure 1.
a_grid=linspace(Params.alowerbar,Params.aupperbar,n_a)';

% Set grid for interest rate, r
r_grid=linspace(-0.05,0.05,n_r)'; % This grid is substantially wider than the actual likely equilibrium values and so is somewhat overkill.
% % Set grid for tax rate
% tau_grid=linspace(); % Can calculate this from the gov budget constraint

%Bring model into the notational conventions used by the toolkit
d_grid=linspace(0,1,n_d)'; % Labor supply
p_grid=r_grid;

n_z=n_theta;
n_p=n_r;

%Create descriptions of SS values as functions of d_grid, a_grid, s_grid &
%pi_s (used to calculate the integral across the SS dist fn of whatever
%functions you define here)
FnsToEvaluateParamNames(1).Names={};
FnsToEvaluateFn_1 = @(d_val, aprime_val,a_val,z_val) a_val; % Aggregate assets (which is this periods state)
%FnsToEvaluateParamNames(2).Names={'v'};
%FnsToEvaluateFn_2 = @(d_val, aprime_val,a_val,z_val,v) v*(z_val==0); % Total unemployment benefits
FnsToEvaluate={FnsToEvaluateFn_1}; %, FnsToEvaluateFn_2};

%Now define the functions for the General Equilibrium conditions
    %Should be written as LHS of general eqm eqn minus RHS, so that 
    %the closer the value given by the function is to zero, the closer 
    %the general eqm condition is to holding.
GeneralEqmEqnParamNames(1).Names={'B'};
GeneralEqmEqn_1 = @(AggVars,p,B) AggVars(1)-B; %The requirement that the aggregate assets (lending and borrowing) equal zero
% GeneralEqmEqnParamNames(2).Names={'B'};
% GeneralEqmEqn_2 = @(AggVars,p,B) p(2)-AggVars(2)-(p(1)/(1+p(1)))*B; % Government budget constraint
GeneralEqmEqns={GeneralEqmEqn_1}; %, GeneralEqmEqn_2};


%% 
DiscountFactorParamNames={'beta'};

ReturnFn=@(d_val, aprime_val, a_val, z_val,r, gamma, psi, eta, phi, v) GuerrieriLorenzoni2017_ReturnFn(d_val, aprime_val, a_val, z_val,r, gamma, psi, eta, phi, v);
ReturnFnParamNames={'r', 'gamma', 'psi', 'eta', 'phi', 'v'}; %It is important that these are in same order as they appear in 'Aiyagari1994_ReturnFn'

%% Solve

V0=ones(n_a,n_z,'gpuArray');
%Use the toolkit to find the equilibrium price index
GEPriceParamNames={'r'}; %,'tau'

disp('Calculating price vector corresponding to the stationary eqm')
% tic;
heteroagentoptions.verbose=1;
heteroagentoptions.pgrid=p_grid;
[p_eqm,p_eqm_index, MarketClearance]=HeteroAgentStationaryEqm_Case1(V0, n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
% findeqmtime=toc
% Params.q=p_eqm;
save ./SavedOutput/GuerrieriLorenzoni2017.mat p_eqm p_eqm_index MarketClearance

%% %Now that we know what the equilibrium price is, lets calculate a bunch of
load ./SavedOutput/GuerrieriLorenzoni2017.mat p_eqm p_eqm_index MarketClearance

% %other things associated with the equilibrium
p_eqm_index
disp('Calculating various equilibrium objects')
[~,Policy]=ValueFnIter_Case1(V0, n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames);

% PolicyValues=PolicyInd2Val_Case1(Policy,n_d,n_a,n_s,d_grid,a_grid, Parallel);

StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z, simoptions);

AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate,Params, FnsToEvaluateParamNames,n_d, n_a, n_z, d_grid, a_grid,z_grid,pi_z, Parallel);

eqm_MC=real(GeneralEqmConditions_Case1(AggVars,Params.q, GeneralEqmEqns, Params, GeneralEqmEqnParamNames));

save ./SavedOutput/GuerrieriLorenzoni2017.mat p_eqm p_eqm_index MarketClearance Policy StationaryDist


%% Display some output about the solution
GraphString=['Huggett1993_MarketClearance_mu', num2str(Params.mu), 'alowerbar', num2str(abs(Params.alowerbar))];
GraphString=strrep(GraphString, '.', '_');
% Create a graph displaying the market clearance
fig1=figure(1);
plot(p_grid,MarketClearance, 'x',p_grid, zeros(n_p,1),p_eqm, 0, 'r o')
title(['Market clearance: mu=', num2str(Params.mu), ' alowerbar=', num2str(Params.alowerbar)],'FontSize',18);
xlabel('p','FontSize',16); ylabel('tilde(p)-p','FontSize',16);
set(fig1, 'Color', 'white');     % white bckgr
set(fig1, 'Unit', 'centimeters');  % set unit to cm
set(fig1,'position',[0 0 20 10]);  % set size
export_fig(fig1, ...            % figure handle
    ['./SavedOutput/Graphs/',GraphString,'.pdf'],... % name of output file without extension   %  '-painters', ...            % renderer
    '-pdf', ...                 % file format
    '-r300' );                  % resolution

%plot(cumsum(sum(StationaryDist,2))) %Plot the asset cdf

fprintf('For parameter values mu=%.2f, alowerbar=%.2f \n', [Params.mu,Params.alowerbar])
fprintf('The table elements are q=%.4f, r=%.4f \n',[Params.q,100*Params.r])

%% Outputs of the function
% Tables: q,r

Output.q=gather(Params.q);
Output.r=gather(100*Params.r);
Output.Policy=gather(Policy);
Output.StationaryDist=gather(StationaryDist);
Output.a_grid=gather(a_grid);

end




