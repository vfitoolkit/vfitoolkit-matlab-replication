function F=ReturnFn_HIP(aprime,a,betahat,zhat,v,workret_j,exper_j,Wbar_j,Pb,crra,Phibar,Ybar_T,exper_T,alpha_i)
% Guvenen (2007) HIP model with Bayesian learning, KNOWN-ALPHA version (NOTES.md stage c).
% Exogenous states: predictive beliefs (betahat, zhat) and the income innovation v.
% alpha_i is the (known) intercept of the current alpha-type, passed as a parameter.
% Realized log income is y_t = alpha_i + betahat*exper_j + zhat + v, so
%   Working: Y = exp(alpha_i + betahat*exper_j + zhat + v)
%   Retired: pension = Phibar * PhiFn(Ytilde_T) * Ybar_T with Ytilde_T = Y_T/Ybar_T and
%     Y_T = exp(alpha_i + betahat*exper_T + zhat + v) from the FROZEN last-working-period
%     states (their transitions are identity from the last working period on).
% Age-dependent parameters: workret_j, exper_j, Wbar_j (limit on this period's aprime).
% Budget: c + Pb*aprime = a + Y;  constraint aprime >= Wbar_j.

F=-Inf;

if workret_j==1
    Y=exp(alpha_i+betahat*exper_j+zhat+v);
else
    Ytilde=exp(alpha_i+betahat*exper_T+zhat+v)/Ybar_T;
    if Ytilde<0.3
        PhiFn=0.9*Ytilde;
    elseif Ytilde<=2
        PhiFn=0.27+0.32*(Ytilde-0.3);
    elseif Ytilde<=4.1
        PhiFn=0.81+0.15*(Ytilde-2);
    else
        PhiFn=1.1;
    end
    Y=Phibar*PhiFn*Ybar_T;
end

c=a+Y-Pb*aprime;

if c>0 && aprime>=Wbar_j-1e-9
    F=(c^(1-crra))/(1-crra);
end

end
