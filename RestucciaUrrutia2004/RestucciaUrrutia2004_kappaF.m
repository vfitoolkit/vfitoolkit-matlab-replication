function kappaF=RestucciaUrrutia2004_kappaF(bhatprime,sprime,h,bhat,s,b,thetahat,earningshock,agej,w,kappa0,kappa1,nlowerbar,nupperbar,f,psi0,psi1)
% Case 2, Finite Horizon.
% bhatprime: (decision variable), acquired ability of child (only j==1)
% h: (endogenous state), human capital (both j==1 & 2)
% bhat: (endogenous state), acquired ability of child (only j==2)
% theta: (endogenous state), completed or only partial college (only j==2)
% b: (exogenous state) innate ability of child (both j==1 & 2)
% Note: I use theta to take different values to those used by Restuccia & Urrutia (2005)
% Note: The VFI toolkit is used to shrink the 'grid' of these variables to just a single point for the ages in which they are irrelevant.

% agej: age (takes values of 1,2)

kappaF=0;
if agej==2
    q_bhat=min(psi0*(1+bhat)^psi1,1);
    if s==1 && thetahat>q_bhat % Partial college
        Kappa_wh2=max(kappa1-kappa0*w*h,0);
        kappaF=Kappa_wh2*f*nlowerbar;
    elseif s==1 && thetahat<=q_bhat % Completed college
        Kappa_wh2=max(kappa1-kappa0*w*h,0);
        kappaF=Kappa_wh2*f*nupperbar;
    end
end

end