function c=RestucciaUrrutia2004_Consumption(bhatprime,sprime,h,bhat,s,b,thetahat,earningsshock,j,w,g,gamma,tau,kappa0,kappa1,plowerbar, pupperbar,nlowerbar,nupperbar,f,psi0,psi1)
% Case 2, Finite Horizon.
% bhatprime: (decision variable), acquired ability of child (only j==1)
% h: (endogenous state), human capital (both j==1 & 2)
% bhat: (endogenous state), acquired ability of child (only j==2)
% theta: (endogenous state), completed or only partial college (only j==2)
% b: (exogenous state) innate ability of child (both j==1 & 2)
% Note: I use theta to take different values to those used by Restuccia & Urrutia (2005)
% Note: The VFI toolkit is used to shrink the 'grid' of these variables to just a single point for the ages in which they are irrelevant.

% j: age (takes values of 1,2)

% Return function of young.
c=0; % Just a placeholder
if j==1
    % bhat=G(b,e,g) and G(b,e,g)=b*(e+g)^gamma
    % Rearrange these to get e,
    % (Done in this way since we want bhatprime to stay on the
    % bhat grid and this is the easiest way to ensure this)
    e=(bhatprime/b)^(1/gamma)-g;
    c=(1-tau)*w*h*earningsshock-e;
elseif j==2
    q_bhat=min(psi0*(1+bhat)^psi1,1);
    if s==0 % No college
        c=(1-tau)*(w*h*earningsshock+w*bhat); % Note: h1'=bhat (in notation of Restuccia & Urrutia, 2005)
    elseif s==1 && thetahat>q_bhat % Partial college
        Kappa_wh2=max(kappa1-kappa0*w*h*earningsshock,0);
        c=(1-tau)*(w*h*earningsshock+w*(plowerbar*bhat)*(1-nlowerbar))-(1-Kappa_wh2)*f*nlowerbar; % Note: h1'=plowerbar*bhat (in notation of Restuccia & Urrutia, 2005)
    elseif s==1 && thetahat<=q_bhat % Completed college
        Kappa_wh2=max(kappa1-kappa0*w*h*earningsshock,0);
        c=(1-tau)*(w*h*earningsshock+w*(pupperbar*bhat)*(1-nupperbar))-(1-Kappa_wh2)*f*nupperbar; % Note: h1'=pupperbar*bhat (in notation of Restuccia & Urrutia, 2005)
    end
end

end