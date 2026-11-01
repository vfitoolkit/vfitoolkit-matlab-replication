function F=AV2001_HouseholdReturnFn(aprime,a,z,i,w,alpha,tau,eta)
% Household problem with an EXOGENOUS search intensity eta (Tables 2 and 4).
% Eqs (3.3), (3.5), (3.6) of AV2001, with the unemployed budget from eq (13)
% of the 1998 WP (no wage income while unemployed).
%
% z encodes both the agent's status and their income this period (in units of w):
%    z = 10*status + inc
%    status: 0 employed (inc=1), 1 unemployed (inc=lambda in the period of
%    firing, otherwise 0), 2 retired (inc=0), 3 dead (absorbing, value zero).
% Only the transition matrix knows which establishment (e,s) an employed agent
% is at; the return function does not need it.
%
% Dead state: return is zero for every aprime, so V(a,dead)=0 and retirees get
% beta*[sigma*R(a')+(1-sigma)*0] with no bequest motive. Assets are inherited
% by the offspring, which the stationary distribution handles with a different
% transition matrix (retired -> unemployed) than the value function uses.

status=floor(z/10);
inc=z-10*status;

c=(1+i)*a+w*inc-aprime;

F=-Inf;
if c>0
    F=log(c);
    if status==1
        F=F+alpha*((1-eta)^tau-1)/tau; % u(1-eta), eq (2.2)
    end
end
if status==3
    F=0;
end

end
