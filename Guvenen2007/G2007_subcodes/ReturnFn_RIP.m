function F=ReturnFn_RIP(aprime,a,alpha,z,e,workret_j,exper_j,Wbar_j,betabar,Pb,crra,Phibar,Ybar_T,exper_T)
% Guvenen (2007) RIP model return function.
% States: assets a; alpha (fixed effect, identity transition); z (AR(1), frozen in
% retirement); e (transitory, iid while working, frozen in retirement).
% Age-dependent parameters: workret_j (1=working), exper_j (experience, 0 at age 25),
% Wbar_j (natural borrowing limit applying to the aprime chosen this period).
%
% Working: Y = exp(alpha + betabar*exper_j + z + e)
% Retired: pension = Phibar * PhiFn(Ytilde_T) * Ybar_T, where Ytilde_T = Y_T/Ybar_T,
%   Y_T = exp(alpha + betabar*exper_T + z + e) using the FROZEN (z,e) from the last
%   working period (their transitions are identity from the last working period on),
%   and PhiFn is the Storesletten-Telmer-Yaron piecewise-linear replacement rule
%   (Guvenen 2007, p.701), Phibar ~ 1/1.40 (footnote 18).
%
% Budget: c + Pb*aprime = a + Y;  constraint aprime >= Wbar_j.

F=-Inf;

if workret_j==1
    Y=exp(alpha+betabar*exper_j+z+e);
else
    Ytilde=exp(alpha+betabar*exper_T+z+e)/Ybar_T;
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
