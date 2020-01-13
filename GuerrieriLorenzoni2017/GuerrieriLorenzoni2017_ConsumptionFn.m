function c=GuerrieriLorenzoni2017_ConsumptionFn(d_val, aprime_val, a_val, z_val,r, v, B, Bprime)

% Guerrieri & Lorenzoni (2017) refer to: d as n_it, aprime_val as bit+1, a_val as b_it, and z_val as theta_it.

% In steady-state we know that tau=u*v+r/(1+r)*B
tau=0.0607*v+B-(1/(1+r))*Bprime; % From pg 1434, notice that q_t=1/(1+r)
% In steady-state we know that tau=u*v+r/(1+r)*B, as Bprime=B, so the formula could be simplified.
% But this more general form will allow us to later do the 'Fiscal Policy' transitions.
% Tax and transfers:
tau_tilde=tau-v*(z_val==0); % Note: by definition z cannot actually be negative

c=a_val+z_val*d_val-tau_tilde-(1/(1+r))*aprime_val; % q=1/(1+r)
% if c>0
%     % Note: this won't work if gamma=1 or eta=1
%     F=(c^(1-gamma) -1)/(1-gamma)+psi*(c^(1-eta) -1)/(1-eta);
% end
% 
% % Impose borrowing constraint
% if aprime_val<-phi 
%     F=-Inf;
% end

end