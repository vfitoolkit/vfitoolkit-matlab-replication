function F=GuerrieriLorenzoni2017_StickyWageReturnFn(d_val, aprime_val, a_val, z_val,r, gamma, psi, eta, phi, v, omega)%, tau)

% Guerrieri & Lorenzoni (2017) refer to: d as n_it, aprime_val as bit+1, a_val as b_it, and z_val as theta_it.

% In steady-state we know that tau=u*v+r/(1+r)*B
tau=0.0607*0.1+(r/(1+r))*1.6;

% Tax and transfers:
tau_tilde=tau-v*(z_val==0); % Note: by definition z cannot actually be negative

F=-Inf;
c=a_val+z_val*d_val-tau_tilde-(1/(1+r))*aprime_val; % q=1/(1+r)
if c>0
    % Note: this won't work if gamma=1 or eta=1
    F=(c^(1-gamma) -1)/(1-gamma)+(1/(1-omega))*psi*((1-d_val)^(1-eta) -1)/(1-eta);
end

% Impose borrowing constraint
if aprime_val<-phi 
    F=-Inf;
end

end