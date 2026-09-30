function k=Kitao2008_kFn(aprime,eprime,a,e,eta,theta,d,upsilon1,upsilon2,r,w,delta,phi,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k)
% Capital used by entrepreneur. Returns 0 if worker.

k=0; % just to make GPU happy

if e==1 % Entrepreneurs' problem (for workers, k stays at zero)
    % The entrepreneurs production problem is static
    [k,n_unused,y_unused]=Kitao2008_StaticEntrepreneurProblem(a,eta,r,w,d,phi,delta,theta,upsilon1,upsilon2,tau_a0,tau_a1,tau_a2,tau_I,taxincome,tau_k);
end

end
