function c=ImrohorogluImrohorogluJoines1995_ConsumptionFn(kprime,k,z,r,tau_u, tau_s,h,zeta,epsilon_j,I_j,alpha,delta, A,SSdivw, Tr_beq) %i

w=(1-alpha)*(A^(1/(1-alpha)))*((r+delta)/alpha)^(alpha/(alpha-1));
% Yhat=(rhat+delta)/(1-alpha);

earnings=z*w*h*epsilon_j*I_j;
u=zeta*w*h; % Unemployment benefits are at replacement rate phi;
unemploymentbenefit=(1-z)*u*I_j;
socialsecuritybenefit=w*SSdivw*(1-I_j);
% q=earnings+unemploymentbenefit+socialsecuritybenefit; % This is the notation used by IIJ1995 and is just here to ease comprehension. 

c=(1+r)*k+(1-tau_u-tau_s)*earnings+unemploymentbenefit+socialsecuritybenefit+Tr_beq-kprime;

% 
% F=-Inf;
% 
% if c>0
%     F=(c.^(1-gamma))/(1-gamma);
% end
% 
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