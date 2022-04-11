function F=GuerrieriLorenzoni2017_DurableGoods_ReturnFn(d_val, a1prime_val, a2prime_val, a1_val, a2_val, z_val,r, alpha, gamma, psi, eta, phi_k, delta, zeta, chi, v, B, Bprime, omega)

% Guerrieri & Lorenzoni (2017) refer to: d as n_it, a1prime_val as b_it+1, 
% a2prime_val as k_it+1, a1_val as b_it, a2_val as k_it, and z_val as theta_it.
% omega is the 'wedge' induced by sticky wages in the New Keynesian version
% of model in Section 4 of GL (2017), pg 1450. Notice that when omega=0 it
% plays no role and this is the case in the rest of the paper.
% So a2 is durable consumption good, a1 is non-durable consumption good

tau=0.0607*v+B-(1/(1+r))*Bprime; % From pg 1434, notice that q_t=1/(1+r)
% In steady-state we know that tau=u*v+r/(1+r)*B, as Bprime=B, so the formula could be simplified.
% But this more general form will allow us to later do the 'Fiscal Policy' transitions.

% Tax and transfers:
tau_tilde=tau-v*(z_val==0); % Note: by definition z cannot actually be negative

% Adjustment costs g(k_t+1,kt) [pg 48 of online appendix to GL2017]
gk=(a2prime_val-a2_val+delta*a2_val)*(a2prime_val>=a2_val)+((1-zeta)*(a2prime_val-a2_val)+delta*a2_val)*(a2prime_val<a2_val);
% Note: GL2017 choose this specification as it keeps the household problem
% concave, which makes it possible for them to continue to use the kinds of
% numerical solution methods they use (concave problems are easier and
% faster to solve as FOCs are necessary and sufficient). Given how VFI Toolkit works this provides no advantage
% here and could just as well use any other form.

F=-Inf;
% Bottom of page 48 of GL2017 Appendix A.3, we get the following budgent constraint
c=a1_val+z_val*d_val-tau_tilde-(1-chi)*(a1prime_val<0)*(1/(1+r))*a1prime_val-(a1prime_val<0)*(1/(1+r))*a1prime_val-gk; % q=1/(1+r)
% Note: as well as adding capital adjustment cost gk, also introduce intermediation cost between borrowing and lending cost (1-zeta)*(a1prime_val<0)
% Note: when a1prime_val=0 we just get zero in the above, no need to treat it seperate from the < and >
if c>0 && a2_val>0
    % Note: this won't work if gamma=1 or eta=1
    F=(((c^alpha)*(a2_val^(1-alpha)))^(1-gamma))/(1-gamma)+(1/(1-omega))*psi*((1-d_val)^(1-eta))/(1-eta);
end

% Impose borrowing constraint (now a collateral constraint)
if a1prime_val<(-phi_k*a2prime_val)
    F=-Inf;
end

end