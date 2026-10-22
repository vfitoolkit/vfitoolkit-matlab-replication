function F=HVY2011_Table4_Objective(guess,hip,VARLOGE,COVARLOGE,LAG)
% Equally weighted distance between the variances and autocovariances of log earnings by age in a sample (VARLOGE, COVARLOGE)
% and those implied by the RIP (hip=0) or HIP (hip=1) statistical model of earnings, HVY2011 Table 4.
% Follows the authors' code (data/hvy-programfiles/gmm-estimation.f and gmm-estimation-hip.f, function funk): ages are indexed
% H=1,...,38 (ages 23-60), the moments start at H=4 (age 26), where the persistent component starts (z=eta), and the growth
% rate effect beta_i multiplies H-3. The parameters are in logs, rho=1-0.1*exp(guess(3)) is bounded below by 0.91 (RIP) or 0.40
% (HIP), and cov(alpha,beta)=exp(guess(6))-1.
%   log e_iH = alpha_i + beta_i*(H-3) + z_iH + eps_iH,  z_iH = rho*z_iH-1 + eta_iH,  z_i3=0

RETX=38;
varinn=exp(guess(1));   % variance of the persistent shock, sigma_eta^2
varfixed=exp(guess(2)); % variance of the fixed effect, sigma_alpha^2
corr=1-0.1*exp(guess(3)); % rho
vartemp=exp(guess(4));  % variance of the temporary shock, sigma_epsilon^2
if hip==0
    corr=max(0.91,corr);
    varbeta=0;
    covalphabeta=0;
else
    corr=max(0.40,corr);
    varbeta=exp(guess(5));
    covalphabeta=exp(guess(6))-1;
end

VARP=zeros(RETX,1); % variance of the persistent component
VARLOGE2=zeros(RETX,1);
VARP(4)=varinn;
VARLOGE2(4)=varinn+vartemp+varfixed+varbeta+2*covalphabeta;
for H=5:RETX
    VARP(H)=(corr^2)*VARP(H-1)+varinn;
    VARLOGE2(H)=vartemp+varfixed+VARP(H)+varbeta*((H-3)^2)+2*covalphabeta*(H-3);
end
F=sum((VARLOGE2(4:RETX)-VARLOGE(4:RETX)).^2);
for H=4+LAG:RETX
    for IK=1:LAG
        COVARLOGE2=varfixed+(corr^IK)*VARP(H-IK)+covalphabeta*(2*(H-3)-IK)+varbeta*((H-3)*(H-3-IK));
        F=F+(COVARLOGE2-COVARLOGE(H,IK))^2;
    end
end

end
