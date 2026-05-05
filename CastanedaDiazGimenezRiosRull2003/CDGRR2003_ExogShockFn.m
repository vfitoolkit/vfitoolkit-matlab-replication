function [z_grid,pi_z]=CDGRR2003_ExogShockFn(e2,e3,e4,J,p_eg,p_gg,phi1,phi2,Gamma_ee_12,Gamma_ee_13,Gamma_ee_14,Gamma_ee_21,Gamma_ee_23,Gamma_ee_24,Gamma_ee_31,Gamma_ee_32,Gamma_ee_34,Gamma_ee_41,Gamma_ee_42,Gamma_ee_43)
% Creates the transition matrix on the exogenous shocks.

% We know (conceptually, given what model aims to do) that most of the
% probability mass will likely end up on the diagonals. Hence we
% parameterize this as the 'other elements of the row' [note how the
% missing ones are Gamma_ee_11, Gamma_ee_22, etc.]
% This makes things more stable later on when we come to calibrate the
% model. [Calibrating transition matrices can be tricky, because all the
% elements must be between 0 and 1, and the row must sum to one.]

%% First, build Gamma_ee
Gamma_ee_11=1-Gamma_ee_12-Gamma_ee_13-Gamma_ee_14-p_eg;
Gamma_ee_22=1-Gamma_ee_21-Gamma_ee_23-Gamma_ee_24-p_eg;
Gamma_ee_33=1-Gamma_ee_31-Gamma_ee_32-Gamma_ee_34-p_eg;
Gamma_ee_44=1-Gamma_ee_41-Gamma_ee_42-Gamma_ee_43-p_eg;

Gamma_ee=gpuArray([Gamma_ee_11, Gamma_ee_12, Gamma_ee_13, Gamma_ee_14; Gamma_ee_21, Gamma_ee_22, Gamma_ee_23, Gamma_ee_24; Gamma_ee_31, Gamma_ee_32, Gamma_ee_33, Gamma_ee_34; Gamma_ee_41, Gamma_ee_42, Gamma_ee_43, Gamma_ee_44]);
Gamma_ee=Gamma_ee./(sum(Gamma_ee,2)*ones(1,J,'gpuArray')); %This is a normalization of Gamma_ee into a probability matrix
% Gamma_ee=Gamma_ee./(1-p_eg); %this would be equivalent to the normalization in the line above

% And build gammastar
e_grid=gpuArray([1;e2;e3;e4]); % e(s), but just for first four
% Calculate the gammastar's
[~,~,~,gammastar]=MarkovChainMoments(e_grid,Gamma_ee);

%% Now we calculate all of the points for Gamma_re
% Step 1
p_51_Step1=gammastar(1)+phi1*gammastar(2)+phi1^2*gammastar(3)+phi1^3*gammastar(4);
p_52_Step1=(1-phi1)*(gammastar(2)+phi1*gammastar(3)+phi1^2*gammastar(4));
p_53_Step1=(1-phi1)*(gammastar(3)+phi1*gammastar(4));
p_54_Step1=(1-phi1)*gammastar(4);
p_61_Step1=(1-phi1)*gammastar(1);
p_62_Step1=phi1*gammastar(1)+gammastar(2)+phi1*gammastar(3)+phi1^2*gammastar(4);
p_63_Step1=(1-phi1)*(gammastar(3)+phi1*gammastar(4));
p_64_Step1=(1-phi1)*gammastar(4);
p_71_Step1=(1-phi1)*gammastar(1);
p_72_Step1=(1-phi1)*(phi1*gammastar(1)+gammastar(2));
p_73_Step1=phi1^2*gammastar(1)+phi1*gammastar(2)+gammastar(3)+phi1*gammastar(4);
p_74_Step1=(1-phi1)*gammastar(4);
p_81_Step1=(1-phi1)*gammastar(1);
p_82_Step1=(1-phi1)*(phi1*gammastar(1)+gammastar(2));
p_83_Step1=(1-phi1)*(phi1^2*gammastar(1)+phi1*gammastar(2)+gammastar(3));
p_84_Step1=phi1^3*gammastar(1)+phi1^2*gammastar(2)+phi1*gammastar(3)+gammastar(4);
% Step 2
p_51=p_51_Step1+phi2*p_52_Step1+phi2^2*p_53_Step1+phi2^3*p_54_Step1;
p_61=p_61_Step1+phi2*p_62_Step1+phi2^2*p_63_Step1+phi2^3*p_64_Step1;
p_71=p_71_Step1+phi2*p_72_Step1+phi2^2*p_73_Step1+phi2^3*p_74_Step1;
p_81=p_81_Step1+phi2*p_82_Step1+phi2^2*p_83_Step1+phi2^3*p_84_Step1;
p_52=(1-phi2)*(p_52_Step1+phi2*p_53_Step1+phi2^2*p_54_Step1);
p_62=(1-phi2)*(p_62_Step1+phi2*p_63_Step1+phi2^2*p_64_Step1);
p_72=(1-phi2)*(p_72_Step1+phi2*p_73_Step1+phi2^2*p_74_Step1);
p_82=(1-phi2)*(p_82_Step1+phi2*p_83_Step1+phi2^2*p_84_Step1);
p_53=(1-phi2)*(p_53_Step1+phi2*p_54_Step1);
p_63=(1-phi2)*(p_63_Step1+phi2*p_64_Step1);
p_73=(1-phi2)*(p_73_Step1+phi2*p_74_Step1);
p_83=(1-phi2)*(p_83_Step1+phi2*p_84_Step1);
p_54=(1-phi2)*p_54_Step1;
p_64=(1-phi2)*p_64_Step1;
p_74=(1-phi2)*p_74_Step1;
p_84=(1-phi2)*p_84_Step1;

% Now put these into the matrix
Gamma_re=gpuArray([p_51, p_52, p_53, p_54; p_61, p_62, p_63, p_64; p_71, p_72, p_73, p_74; p_81, p_82, p_83, p_84]);
Gamma_re=Gamma_re./(sum(Gamma_re,2)*ones(1,J,'gpuArray')); %This is a normalization of Gamma_re into a probability matrix

%% Create Gamma (pi_z), and we just set up z_grid as indexes 
% (have to use indexes, rather than use e_grid, so that we can track birth/death and ageing for some model statistics)

z_grid=linspace(1,2*J,2*J)'; % age (& determines retirement)

pi_z=[Gamma_ee.*(1-p_eg), diag(p_eg*ones(J,1,'gpuArray')); Gamma_re.*(1-p_gg), diag(p_gg*ones(J,1,'gpuArray'))];  %transmatix is (z,zprime) %dim N_s-by-N_s (s by sprime)



end