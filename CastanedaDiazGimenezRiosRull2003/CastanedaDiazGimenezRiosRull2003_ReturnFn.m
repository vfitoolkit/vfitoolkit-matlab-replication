function F=CastanedaDiazGimenezRiosRull2003_ReturnFn(d1_val, d2_val, a_val, z_val,r,sigma1,sigma2,chi,elle,theta,delta,e1,e2,e3,e4,omega,a0,a1,a2,a3)
% You can include comments in the return function file, but only on
% seperate lines to the commands. If you try putting a comment on the same
% line as a command (and following line is not blank) then this will cause
% an error with the GPU codes.

w=(1-theta)*(((r+delta)/(theta))^(theta/(theta-1)));


% Determine values of e and omega_e as functions of 'age' (z_val)
e=0;
omega_s=0;
if z_val==1
    e=e1;
elseif z_val==2
    e=e2;
elseif z_val==3
    e=e3;
elseif z_val==4
    e=e4;
else
    omega_s=omega;
end

% Calculate income and then income after tax
y=a_val*r+e*d1_val*w+omega_s; %eqn 5

% Calculate income tax
tau_y=a0*(y-a2)+a3*y;
if a1>0
    tau_y=a0*(y-(y^(-a1)+a2)^(-1/a1))+a3*y;
end

% Estate tax is dealt with as part of 'phiaprime' (Case 2 value function problem).

F=-Inf;
c=a_val+y-tau_y-d2_val; % eqn 4
if c>0
    if sigma1==1
        F=log(c)+chi*((elle-d1_val)^(1-sigma2))/(1-sigma2);
    else
        F=(c^(1-sigma1))/(1-sigma1)+chi*((elle-d1_val)^(1-sigma2))/(1-sigma2);
    end
end

end