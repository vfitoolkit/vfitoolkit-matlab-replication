function F=GuerrieriLorenzoni2017_ReturnFn(d, aprime, a, z,r, gamma, psi, eta, phi, v,B,Bprime,omega)%, tau)

% Guerrieri & Lorenzoni (2017) refer to: d as n_it, aprime as b_it+1, a as b_it, and z as theta_it.
% omega is the 'wedge' induced by stick wages in the New Keynesian version
% of model in Section 4 of GL (2017), pg 1450. Notice that when omega=0 it
% plays no role and this is the case in the rest of the paper.

tau=0.0607*v+B-(1/(1+r))*Bprime; % From pg 1434, notice that q_t=1/(1+r)
% In stationary eqm we know that tau=u*v+r/(1+r)*B, as Bprime=B, so the formula could be simplified.
% But this more general form will allow us to later do the 'Fiscal Policy' transitions.

% (Net) Tax and transfers:
tau_tilde=tau-v*(z==0); % Note: by definition z cannot actually be negative

% Note: Agent can choose labour supply, but only gets unemployment benefits
% if z==0, do not get unemployment benefits if choose zero labour supply (d==0).

F=-Inf;
c=a+z*d-tau_tilde-(1/(1+r))*aprime; % q=1/(1+r)
if c>0
    % Note: this won't work if gamma=1 or eta=1
    F=(c^(1-gamma))/(1-gamma)+(1/(1-omega))*psi*((1-d)^(1-eta))/(1-eta);
end

% Impose borrowing constraint
if aprime<-phi 
    F=-Inf;
end

end