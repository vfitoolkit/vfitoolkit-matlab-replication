function income=ImrohorogluImrohorogluJoines1995_IncomeFn(kprime,k,z,r,h,zeta, epsilon_j,I_j,alpha,delta, A,SSdivw, Tr_beq,workinglifeincome,g,agej) %i

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

% income=r*k+earnings+unemploymentbenefit+socialsecuritybenefit+Tr_beq;

% From Figure 3 it became clear that above commented out formula was not what
% is being plotted as the concept plotted does not include capital income
% but does include the pension (can be seen from the shape of the income
% profile during retirement)    
income=earnings+unemploymentbenefit+socialsecuritybenefit+Tr_beq;

end