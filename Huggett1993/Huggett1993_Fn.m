function Output=Huggett1993_Fn(n_a,n_e,n_q,Params,Parallel, vfoptions, simoptions)
%This code replicates the results of Huggett (1993) - Uninsured Idiosyncratic Risk and Aggregate Saving
disp('Running Huggett1993_Fn')
%% Set up

%Create markov process for the exogenous income (based on idea of employment and unemployment states, following Imrohoroglu, 1989).
e_grid=[Params.el,Params.eh]';
pi_e=[1-Params.piehel,Params.piehel; 1-Params.pieheh, Params.pieheh];


%% Grids
% Set grid for asset holdings
Params.aupperbar=4; % Not clear exactly what value Huggett used, have gone with this based on axes of Figure 1.
a_grid=linspace(Params.alowerbar,Params.aupperbar,n_a)';

%Set grid for asset prices, q
q_grid=sort(linspace(0.9,1.1,n_q))'; % This grid is substantially wider than the actual likely equilibrium values and so is somewhat overkill.

%Bring model into the notational conventions used by the toolkit
d_grid=0; %There is no d variable
z_grid=e_grid;
pi_z=pi_e;
p_grid=q_grid;

n_d=0;
n_z=n_e;
n_p=n_q;

%Create descriptions of SS values as functions of d_grid, a_grid, s_grid &
%pi_s (used to calculate the integral across the SS dist fn of whatever
%functions you define here)
SSvalueParamNames(1).Names={};
SSvaluesFn_1 = @(aprime_val,a_val,s_val) a_val; %We just want the aggregate assets (which is this periods state)
SSvaluesFn={SSvaluesFn_1};

%Now define the functions for the General Equilibrium conditions
    %Should be written as LHS of general eqm eqn minus RHS, so that 
    %the closer the value given by the function is to zero, the closer 
    %the general eqm condition is to holding.
GeneralEqmEqnParamNames(1).Names={};
GeneralEqmEqn_1 = @(AggVars,p) AggVars; %The requirement that the aggregate assets (lending and borrowing) equal zero
GeneralEqmEqns={GeneralEqmEqn_1};

disp('sizes')
n_a
n_z
n_p


%% 
DiscountFactorParamNames={'beta'};

ReturnFn=@(aprime_val, a_val, s_val,mu,q) Huggett1993_ReturnFn(aprime_val, a_val, s_val,mu,q);
ReturnFnParamNames={'mu','q'}; %It is important that these are in same order as they appear in 'Aiyagari1994_ReturnFn'

%% Solve

V0=ones(n_a,n_z,'gpuArray'); %(a,s)
%Use the toolkit to find the equilibrium price index
GEPriceParamNames={'q'};

disp('Calculating price vector corresponding to the stationary eqm')
% tic;
heteroagentoptions.pgrid=p_grid;
[p_eqm,p_eqm_index, MarketClearance]=HeteroAgentStationaryEqm_Case1(V0, n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, SSvaluesFn, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, SSvalueParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);
% findeqmtime=toc
Params.q=p_eqm;
save ./SavedOutput/Huggett1993Market.mat p_eqm p_eqm_index MarketClearance

% Equilibrium interest rate (at annual rate; model period is 1/6th of a year)
Params.r=(1+(1-Params.q)/Params.q)^6-1;

% %Now that we know what the equilibrium price is, lets calculate a bunch of
% %other things associated with the equilibrium
p_eqm_index
disp('Calculating various equilibrium objects')
[~,Policy]=ValueFnIter_Case1(V0, n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames);

% PolicyValues=PolicyInd2Val_Case1(Policy,n_d,n_a,n_s,d_grid,a_grid, Parallel);

StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z, simoptions);

SSvalues_AggVars=SSvalues_AggVars_Case1(StationaryDist, Policy, SSvaluesFn,Params, SSvalueParamNames,n_d, n_a, n_z, d_grid, a_grid,z_grid,pi_z, Parallel);

eqm_MC=real(GeneralEqmConditions_Case1(SSvalues_AggVars,Params.q, GeneralEqmEqns, Params, GeneralEqmEqnParamNames));

save ./SavedOutput/Huggett1993SSObjects.mat p_eqm Policy StationaryDist


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




