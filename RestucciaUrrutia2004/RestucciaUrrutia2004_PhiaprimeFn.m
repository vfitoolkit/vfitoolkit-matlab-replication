function Phi=RestucciaUrrutia2004_PhiaprimeFn(bhatprime,sprime,h,bhat,s,b,thetahat,earningsshock,agej,plowerbar,pupperbar,zeta,psi0,psi1,log_lb_h,log_ub_h,loghgridspacing,n_h,log_lb_bhat,log_ub_bhat,logbhatgridspacing,n_bhat)
% Case 2, Finite Horizon. (d,a,z) % Note that j=1 requires d, and j=2
% requires (a,z). But just has to use (d,a,z) as allowing this aspect to
% depend on age is to complex for the moment. (Conceptually easy, but would
% require some thought on how to code.)
% bhatprime: (decision variable), acquired ability of child (only j==1)
% sprime: (decision variable), send child to university (only j==1)
% h: (endogenous state), human capital (both j==1 & 2)
% bhat: (endogenous state), acquired ability of child (only j==2)
% s: (endogenous state), went to high school or university (only j==2)
% b: (exogenous state) innate ability of child (both j==1 & 2)
% thetahat: (exogenous state), (determines whether) completed or only partial college (only j==2)
% earningsshock: (exogenous state), only used for Table 7
% Note: The VFI toolkit is used to shrink the 'grid' of these variables to just a single point for the ages in which they are irrelevant.

% j: age (takes values of 1,2)

% The 'gridspacing' variables (loghgridspacing and logbhatgridspacing) are both
% already measured in terms of the log grid, hence they do not need log
% here.

Phivec1=0; Phivec2=0; Phivec3=0; %(hprime,bhatprime,sprime)
hprime=0;
if agej==1
    % Next period (j==2) value of h is given by
    loghprime=log(zeta*h); %(note that h_grid is equally spaced in logs)
    % Find index of nearest point to hprime on h_grid
    if loghprime<=log_lb_h % If goes off bottom of grid
        Phivec1=1;
    elseif loghprime>=log_ub_h % If goes off top of grid
        Phivec1=n_h;
    else
        Phivec1=1+round((loghprime-log_lb_h)/loghgridspacing);
%         htemp=(hprime-log_lb_h)/loghgridspacing;
%         h_rem=htemp-floor(htemp);
%         if h_rem<=loghgridspacing/2 % Note that this trick requires h_grid to be both uniformly spaced and strictly positive
%             Phivec1=1+floor(htemp);
%         else
%             Phivec1=1+ceil(htemp);
%         end
    end
    % Next period (j==2) value of bhat is given by
    % Find index of nearest point to bhatprime on bhat_grid
    % Since bhat_grid is equally spaced in logs this is done in terms of logs.
    logbhatprime=log(bhatprime);
    if logbhatprime<=log_lb_bhat % If goes off bottom of grid
        Phivec2=1;
    elseif logbhatprime>=log_ub_bhat % If goes off top of grid
        Phivec2=n_bhat;
    else
%         bhattemp=(logbhatprime-log_lb_bhat)/logbhatgridspacing;
        Phivec2=1+round((logbhatprime-log_lb_bhat)/logbhatgridspacing);
%         bhat_rem=bhattemp-floor(bhattemp);
%         if bhat_rem<=logbhatgridspacing/2 % Note that this trick requires h_grid to be both uniformly spaced and strictly positive
%             Phivec2=1+floor(bhattemp);
%         else
%             Phivec2=1+ceil(bhattemp);
%         end
    end
    % Next period (j==2) value of s is given by
    Phivec3=sprime+1; % This one is trivial (sprime takes values 0 and 1)
else %if agej==2
    % Next period (j==1) value of h is given by
    q_bhat=min(psi0*(1+bhat)^psi1,1);
    if s==0
        hprime=bhat; 
    elseif s==1 && thetahat>q_bhat
        hprime=plowerbar*bhat;
    elseif s==1 && thetahat<=q_bhat 
        hprime=pupperbar*bhat;
    end
    loghprime=log(hprime); %(note that h_grid is equally spaced in logs)
    % Find index of nearest point to hprime on h_grid
    if loghprime<=log_lb_h % If goes off bottom of grid
        Phivec1=1;
    elseif loghprime>=log_ub_h % If goes off top of grid
        Phivec1=n_h;
    else
        Phivec1=1+round((loghprime-log_lb_h)/loghgridspacing);
%         htemp=abs((hprime-log_lb_h))/loghgridspacing;
%         h_rem=htemp-floor(htemp);
%         if h_rem<=loghgridspacing/2 % Note that this trick requires (log of) h_grid to be both uniformly spaced and strictly positive
%             Phivec1=1+floor(htemp);
%         else
%             Phivec1=1+ceil(htemp);
%         end
    end
    % Next period (j==1) value of bhat is given by
    Phivec2=1; % Variable is not relevant to age j==1, and so is 'made to disappear' by setting it to just a single point.    
    % Next period (j==1) value of theta is given by
    Phivec3=1; % Variable is not relevant to age j==1, and so is 'made to disappear' by setting it to just a single point.
end

% jplusone_ind=rem(agej,2)+1;
% Phi=Phivec1+(Phivec2-1)*n_h*jplusone_ind+(Phivec3-1)*n_h*jplusone_ind*n_bhat*jplusone_ind;
% Is just implementing sub2ind_homemade([n_h,n_bhat,n_theta],Phivec)

Phi=Phivec1+(Phivec2-1)*n_h+(Phivec3-1)*n_h*n_bhat;
% Is just implementing sub2ind_homemade([n_h,n_bhat,n_theta],Phivec)

end