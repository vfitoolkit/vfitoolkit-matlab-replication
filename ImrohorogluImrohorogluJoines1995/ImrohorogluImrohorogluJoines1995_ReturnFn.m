function F=ImrohorogluImrohorogluJoines1995_ReturnFn(kprime,k,z,r,tau_u, tau_s,gamma, h,zeta, epsilon_j,I_j,alpha,delta, A,SSdivw, Tr_beq,MedicalShock,workinglifeincome,g,agej,LumpSum) %i
% Note: In baseline model MedicalShocks=0, workinglifeincome=1, and g=0 and
% hence none of these inputs would be required if we just wanted to solve
% the baseline model. agej would also not be needed.
% LumpSum is zero. Is just needed for calculating the welfare benefits of various reforms.

w=(1-alpha)*(A^(1/(1-alpha)))*((r+delta)/alpha)^(alpha/(alpha-1));
% Yhat=(rhat+delta)/(1-alpha);

earnings=z*w*h*epsilon_j*I_j;
u=zeta*w*h; % Unemployment benefits are at replacement rate phi;
unemploymentbenefit=(1-z)*u*I_j;
socialsecuritybenefit=w*SSdivw*(1-I_j);
if g>0 % So not using baseline model
    socialsecuritybenefit=w*SSdivw*(1-I_j)*workinglifeincome*(1/((1+g)^agej));
    % Note: the actual social security benefit is constant, but the model
    % has been detrended by (1+g)^t and it is for these reason it here
    % appears to decrease with age.
end
% q=earnings+unemploymentbenefit+socialsecuritybenefit; % This is the notation used by IIJ1995 and is just here to ease comprehension. 

% IIJ1995 includes an extension to 'medical shocks' (see pg 110 of IIJ1995)
medicalexpense=0;
if MedicalShock>0
    % In retirement z now becomes the 'cost of illness' (when non-zero valued, the household is ill, when zero valued the household is healthy)
    % The cost of illness is measured as a fraction of employed wage (which is w*h)
    medicalexpense=(1-I_j)*z*w*h; % Note that this will be zero when working age (as I_j will be 1).
end

c=(1+r)*k+(1-tau_u-tau_s)*earnings+unemploymentbenefit+socialsecuritybenefit+Tr_beq-medicalexpense-kprime+LumpSum;

F=-Inf;

if c>0
    F=(c.^(1-gamma))/(1-gamma);
end

% Borrowing constaint is imposed via grid on assets

% if kprime<0 && jj==N
%     % have to die with non-negative assets
%     F=-Inf;
% end
% 
% if kprime<alowerbar*what
%     % Impose the borrowing limit
%     F=-Inf;
% end

    
end