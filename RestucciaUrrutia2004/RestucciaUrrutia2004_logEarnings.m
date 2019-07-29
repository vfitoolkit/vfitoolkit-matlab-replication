function logwH=RestucciaUrrutia2004_logEarnings(bhatprime,sprime,h,bhat,s,b,thetahat,earningshock,agej,plowerbar, pupperbar,nlowerbar,nupperbar,phi0,phi1,w)
% Case 2, Finite Horizon.
% bhatprime: (decision variable), acquired ability of child (only j==1)
% h: (endogenous state), human capital (both j==1 & 2)
% bhat: (endogenous state), acquired ability of child (only j==2)
% theta: (endogenous state), completed or only partial college (only j==2)
% b: (exogenous state) innate ability of child (both j==1 & 2)
% Note: I use theta to take different values to those used by Restuccia & Urrutia (2005)
% Note: The VFI toolkit is used to shrink the 'grid' of these variables to just a single point for the ages in which they are irrelevant.

% agej: age (takes values of 1,2)

% Return function of young.
H=h;
if agej==2
    q_bhat=min(phi0*(1+bhat)^phi1,1);
    if s==1 && thetahat>q_bhat % Partial college
        H=H+(plowerbar*bhat)*(1-nlowerbar); % Note: h1'=plowerbar*bhat (in notation of Restuccia & Urrutia, 2005)
    elseif s==1 && thetahat<=q_bhat % Completed college
        H=H+(pupperbar*bhat)*(1-nupperbar); % Note: h1'=pupperbar*bhat (in notation of Restuccia & Urrutia, 2005)
    end
end
logwH=log(w*h*earningshock);

end