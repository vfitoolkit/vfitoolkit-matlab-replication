function [eta_grid,pi_eta]=Kitao2008_NoEntrepreneurs2_ExogShockFn(e4)

% Process on exogenous shocks
% From Table 5 of CDGRR2003
e1=1/3.15; e2=1; e3=9.78/3.15;
% Create the grids
eta_grid=[e1,e2,e3,e4]';

% From Table 4 (the diagonal elements of Gamma_ee can be considered as residually determined to make each row of Gamma (Gamma_ee) add up to 1 (to 1-p_eg); 
%    in principle it doesn't matter which element in each row is residual, but from perspective of calibration it is easiest to 
%    use the diagonals to avoid contraint that they as transition probabilities they must be >=0 from binding)
Gamma_ee_12=0.0114; Gamma_ee_13=0.0039; Gamma_ee_14=0.0001; % Gamma_ee_11=0.9624;
Gamma_ee_21=0.0307; Gamma_ee_23=0.0037; Gamma_ee_24=0;       % Gamma_ee_22=0.9433;
Gamma_ee_31=0.015;  Gamma_ee_32=0.0043; Gamma_ee_34=0.0002;  % Gamma_ee_33=0.9582;
Gamma_ee_41=0.1066; Gamma_ee_42=0.0049; Gamma_ee_43=0.0611;  % Gamma_ee_44=0.8051;

Gamma_ee_11=1-Gamma_ee_12-Gamma_ee_13-Gamma_ee_14;
Gamma_ee_22=1-Gamma_ee_21-Gamma_ee_23-Gamma_ee_24;
Gamma_ee_33=1-Gamma_ee_31-Gamma_ee_32-Gamma_ee_34;
Gamma_ee_44=1-Gamma_ee_41-Gamma_ee_42-Gamma_ee_43;

Gamma_ee=[Gamma_ee_11, Gamma_ee_12, Gamma_ee_13, Gamma_ee_14;...
    Gamma_ee_21, Gamma_ee_22, Gamma_ee_23, Gamma_ee_24;...
    Gamma_ee_31, Gamma_ee_32, Gamma_ee_33, Gamma_ee_34;...
    Gamma_ee_41, Gamma_ee_42, Gamma_ee_43, Gamma_ee_44];

pi_eta=Gamma_ee; % Set pi_eta to the transition matrix for working ages from CDGRR2003, which they called Gamme_ee

%% Normalize eta_grid so that the unconditional mean of eta is unity
% Kitao (2008), Section 3.1: "The grid of eta is normalized so that the unconditional mean of eta is
% unity". That is done for the benchmark, and economy 'no entrepreneurs 1' inherits it, so it should
% be done here too, otherwise the three economies are not in comparable units.
% This matters even though scaling eta leaves the prices r and w untouched (scaling eta by lambda
% just scales A, Y, C and G by lambda). The reason is the income tax function: tau_a2 is unit
% dependent, Kitao (2008) notes that if income is scaled by lambda then tau_a2 has to be rescaled to
% tau_a2*lambda^(-tau_a1). Since we hold tau_a2 at the value in her Table 1, the units have to match.
% Note also that this has to be done here inside the ExogShockFn rather than once outside it,
% because e4 is being calibrated (to hit the Gini of wealth) and so the mean of eta moves with it.
% Note: pi_eta does not depend on e4, so the stationary distribution below is always the same.
statdist_eta=ones(1,4)/4;
for ii=1:10^5
    statdist_eta_new=statdist_eta*pi_eta;
    if max(abs(statdist_eta_new-statdist_eta))<10^(-14)
        statdist_eta=statdist_eta_new;
        break
    end
    statdist_eta=statdist_eta_new;
end
eta_grid=eta_grid/(statdist_eta*eta_grid); % Normalize so unconditional mean of eta is unity

end