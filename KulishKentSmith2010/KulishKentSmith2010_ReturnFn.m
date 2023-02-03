function F=KulishKentSmith2010_ReturnFn(l,aprime,a,rho,ej,b1,b2,b3,agej,J,r,delta,A,alpha, g_e)
% g_e=0, is needed for Section 4.6 about labor-augementing productivity growth

F=-Inf;

if agej<=J
    % Rearranging that r=MPK-delta gives the following eqn (MPK is marginal product of capital)
    KdivL=((r+delta)/(alpha*A))^(1/(alpha-1));
    % We know w=MPL (MPL is marginal product of labour)
    w=A*(1-alpha)*(KdivL^alpha); % wage rate (per effective labour unit)
    
    % v_jJ=(b1*agej/J)*normcdf(agej,b2*J,b3*J);
    % Cannot use normcdf() together with arrayfun(), so have to do a manual version
    v_jJ = (b1*agej/J)*(1/2)*(1+erf((agej-b2*J)/(sqrt(2)*b3*J)));
    % Note: normcdf() can be done via erf(): second example under https://au.mathworks.com/help/matlab/ref/erf.html

    c=(1+r)*a+(1-l)*ej*w-(1+g_e)*aprime;

    if c>0
        if rho==1
            F=log(c)+v_jJ*log(l);
        else
            F=(c^(1-rho))/(1-rho)+v_jJ*(l^(1-rho))/(1-rho);
        end
    end

else
    % Because survival is part of the discount factor this is not really needed here
    % It does have to be finite else the 0 of discount factor times this
    % will give an error (i.e., cannot set it to -Inf)
    F=0;
end


end