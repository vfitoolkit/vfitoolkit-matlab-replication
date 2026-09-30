function n=Kitao2008_nFn(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k)
% Labor used by entrepreneur (includes the entrepreneurs own labor endowment eta). Returns 0 if worker.

n=0; % just to make GPU happy

if e==1 % Entrepreneurs' problem (for workers, n stays at zero)
    % The entrepreneurs production problem is static
    [k_unused,n,y_unused]=Kitao2008_StaticEntrepreneurProblem(a,eta,r,w,d,phi,delta,theta,upsilon1,upsilon2,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k);
end

end
