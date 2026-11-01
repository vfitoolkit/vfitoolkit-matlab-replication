function F=AV2001_HouseholdReturnFn_endog(eta,aprime,a,z,i,w,alpha,tau)
% Household problem with ENDOGENOUS search intensity eta (Tables 1, 3, 5).
% eta is the decision variable (the toolkit's d), and it drives the
% semi-exogenous transition of z (job-finding probability eta^xi).
%
% z encodes status and income exactly as in AV2001_HouseholdReturnFn:
%    z = 10*status + inc, status: 0 employed, 1 unemployed, 2 retired, 3 dead.
% Only the unemployed may search: employed and retired agents get -Inf for any
% eta>0 (their transition rows do not depend on eta, so this just pins their
% policy to eta=0). The dead state returns 0 at eta=0, so V(a,dead)=0.

status=floor(z/10);
inc=z-10*status;

c=(1+i)*a+w*inc-aprime;

F=-Inf;
if status==3
    if eta==0
        F=0;
    end
elseif status==1
    if c>0
        F=log(c)+alpha*((1-eta)^tau-1)/tau; % u(1-eta), eq (2.2)
    end
elseif eta==0
    if c>0
        F=log(c);
    end
end

end
