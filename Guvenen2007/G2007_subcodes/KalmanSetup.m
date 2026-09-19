function kf = KalmanSetup(Params)
% Precompute the deterministic Kalman-filter objects for the Guvenen (2007 AER)
% HIP learning problem, working periods t=1..T (t=1 is labor-market entry, age 25).
%
% State: S_t=(alpha, beta, z_t)'
% Observation: y_t = alpha + beta*exper_t + z_t + eps_t,  H_t' = (1, exper_t, 1)
% exper_t = Params.exper1 + (t-1); exper1=0 matches Table 2 column (3), see NOTES.md
%
% Transition: S_{t+1} = F S_t + (0,0,eta_{t+1})', F=diag(1,1,rho), Q=diag(0,0,sigma2_eta)
%
% Prior (Guvenen 2007, p.701): agent knows fraction lambda of the beta variance,
%   P_{1|0} = [ sigma2_alpha,               sqrt(1-lambda)*sigma_ab,  0;
%               sqrt(1-lambda)*sigma_ab,    (1-lambda)*sigma2_beta,   0;
%               0,                          0,                        sigma2_eta ]
% (z_0=0 so the prior on z_1=eta_1 has variance sigma2_eta.)
%
% The recursion is deterministic and common across agents, so everything here is
% computed once: for each t,
%   kf.Ppred(:,:,t) = P_{t|t-1},  kf.Pupd(:,:,t) = P_{t|t}
%   kf.Fv(t) = H_t' P_{t|t-1} H_t + R   (innovation variance, agent's measure: v_t ~ iid N(0,Fv_t))
%   kf.K(:,t) = P_{t|t-1} H_t / Fv(t)   (Kalman gain)
% Belief updating for an individual is then
%   Shat_{t|t} = Shat_{t|t-1} + K_t v_t,   Shat_{t+1|t} = F Shat_{t|t}
%
% Params fields used: rho, sigma2_alpha, sigma2_beta, sigma_ab, sigma2_eta,
%   sigma2_eps, lambda, T, exper1

rho=Params.rho;
F=[1,0,0;0,1,0;0,0,rho];
Q=zeros(3,3); Q(3,3)=Params.sigma2_eta;
R=Params.sigma2_eps;
T=Params.T;

P=[Params.sigma2_alpha, sqrt(1-Params.lambda)*Params.sigma_ab, 0; ...
   sqrt(1-Params.lambda)*Params.sigma_ab, (1-Params.lambda)*Params.sigma2_beta, 0; ...
   0, 0, Params.sigma2_eta];

kf.T=T; kf.F=F; kf.Q=Q; kf.R=R;
kf.P110=P;
kf.exper=Params.exper1+(0:T-1);
kf.H=[ones(1,T); kf.exper; ones(1,T)]; % 3-by-T
kf.Ppred=zeros(3,3,T);
kf.Pupd=zeros(3,3,T);
kf.K=zeros(3,T);
kf.Fv=zeros(1,T);

for t=1:T
    H=kf.H(:,t);
    kf.Ppred(:,:,t)=P;
    Fv=H'*P*H+R;
    K=(P*H)/Fv;
    Pu=P-K*(H'*P);
    Pu=(Pu+Pu')/2; % enforce symmetry against roundoff
    kf.Fv(t)=Fv;
    kf.K(:,t)=K;
    kf.Pupd(:,:,t)=Pu;
    P=F*Pu*F'+Q;
end

end
