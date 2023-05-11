function F=GourioMiao2010_ReturnFn(dividend,kprime,k,z,w,delta,alpha_k,alpha_l,capadjconstant,tau_corp,phi,tau_d,tau_cg)
% Whether we set it up so that dividends or equity issuance is the decision
% variable is unimportant, here I use dividends as the decision variable.

% Note: r is not needed anywhere here, it is relevant to the firm via the discount factor.

F=-Inf;

% We can solve a static problem to get the firm labor input
l=(w/(alpha_l*z*(k^alpha_k)))^(1/(alpha_l-1)); % This is just w=Marg. Prod. Labor, but rearranged

% Output
y=z*(k^alpha_k)*(l^alpha_l);

% Profit
profit=y-w*l; % =(1-alpha_l)*y in the general eqm

% Investment
invest=kprime-(1-delta)*k;

% Capital-adjustment costs
capitaladjcost=(capadjconstant/2)*(invest^2)/k; 

% Taxable corporate income
T=profit-delta*k-phi*capitaladjcost;
% -delta*k: investment expensing
% phi is the fraction of capitaladjcost that can be deducted from corporate taxes (=0 in GM2010 baseline)

% Firms financing constraint gives the new equity issuance
s=dividend+invest+capitaladjcost-(profit-tau_corp*T);

% Firms per-period objective
if s>=0 % enforce that 'no share repurchases allowed'
    F=((1-tau_d)/(1-tau_cg))*dividend-s;
end

% Note: dividend payments cannot be negative is enforced by the grid on dividends which has a minimum value of zero


%% GM2010:
% Define XXX=(z*Anorm*k^alpha)^(1/(1-nu)) *(1-nu) 
% Define YYY=delta*tauc*k+(1-delta)*k-kprime-capitaladjcost -fc*notchangingcap + slowerbar 
% where Anorm=1
%       slowerbar=0, is the repurchase constraint, s<slowerbar
%       nu=0.65, is the exponent on labor in production function (what I call alpha_l)
%
% divifszero = (1-tauc)*XXX*(nu/wage)^(nu/(1-nu)) + YYY;
%    is the dividend if s>=slowerbar; it is also minus the s required if the dividend is zero.
%
% They then define 

% Note: XXX*(nu/wage)^(nu/(1-nu)) is y*(1-nu) which is profits (assumes general eqm wage).
% Hence divifszero is essentially my line 30, impose s=0 and then you get
% dividend=(1-tauc)*profit+((1-delta)*k-kprime)+capitaladjcost-tau_corp(-delta*k-phi*capitaladjcost),
% which is now like their (1-tauc)*XXX*stuff+YYY

% Pretty sure they mess this up as it is only true that you can do the
% profits=(1-nu)*y trick if wage is the optimal wage. But it is typically
% not the optimal wage.
% However, in general eqm it should not matter as long as they have the correct w.

% They find a general eqm wage of 1.2606

end
