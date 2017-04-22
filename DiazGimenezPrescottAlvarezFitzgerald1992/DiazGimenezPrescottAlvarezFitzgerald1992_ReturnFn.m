function F=DiazGimenezPrescottAlvarezFitzgerald1992_ReturnFn(N_val,Aprime_val,Kprime_val, A_val,K_val, s_val,z_val,phi,e,omega,w1,w2,w3,w4,theta,i_d,i_l,mew,alpha,alpha_k,gamma,tau,psi,delta_r,Experiment)
% A is called a by Diaz-Gimenez, Prescott, Alvarez & Fitzgerald (1992) [renamed so can use VFI toolkit standard notation of 'a' being endogenous state vector]

F=-Inf;
C=nan;

%firstly calculate the values of x_d, x_s (note that this implicitly
%imposes constraint given by eqn 23)
x_d=0;
x_s=0;
if K_val>Kprime_val
%     x_d=0;
    x_s=phi*(K_val-Kprime_val);
elseif K_val<Kprime_val
    x_d=Kprime_val-K_val;
%     x_s=0;
end


%secondly calculate the values of d & l (note that this implicitly imposes
%constraint given by eqn 22; paper mixes up e(z) and epsilon(z) here; [e(z)=epsilon(z)-1; both are inflation, one expresses 5% as 1.05 other as 0.05])
%if Aprime_val==0;
d=0;
l=0;
if Aprime_val>0
    d=(1+e)*Aprime_val;
%     l=0;
elseif Aprime_val<0
%     d=0;
    l=-Aprime_val*(1+e); %the minus is because Aprime_val itself is negative
end

pension=0;
if A_val==0 && K_val==0 && s_val==3 %retired, no assets, no house
    pension=omega;
end

w=nan;
if s_val==1
    w=w1;
elseif s_val==2
    w=w2;
elseif s_val==3
    w=w3;
elseif s_val==4
    w=w4;
end
if Experiment==2
    if z_val==1
        w=w*1.02;
    elseif z_val
        w=w*0.98;
    end
end

C=A_val+(1-theta)*(w*N_val+d*i_d-l*i_l) + x_s-x_d +l-d+pension-mew*Kprime_val; % Budget constraint: Eqn 20

if C>0 %If this is not met, then just leave Fmatrix as -Inf
    if s_val==1 || s_val==2
        if K_val==0
            % The following eqn comes from replacing C1 and C2 with C using FOCs from the problem of choosing between c1 and c2 on middle of pg.550
            % Observe that you can rewrite utility fn as (C1^(1-alpha_k/alpha) * C2^(alpha_k/alpha))^alpha
            % As this is just a monotone transformation of the CES utility formulation you get standard result that you spend 1-alpha_k/alpha fraction of C on C1, and alpha_k/alpha fraction of C on C2 at price of
            % gamma per unit. So C1=(1-alpha_k/alpha)C, C2=(1/gamma)*(alpha_k/alpha)*C. Substitute these and rearrange to get the below expression.
            F=( ( C^alpha * (1-alpha_k/alpha)^(alpha-alpha_k) * (alpha_k/(gamma*alpha))^alpha_k * (tau-N_val)^(1-alpha) )^(1-psi) )/(1-psi);
        
            % Impose collateral constraint (eqn 21)
            if l>phi*Kprime_val*(1+e)
                F=-Inf;
            end
        end
        if K_val>0 %ie. if K=3
            F=( ( C^(alpha-alpha_k) * K_val^alpha_k * (tau-N_val)^(1-alpha) )^(1-psi) )/(1-psi);
        end
    elseif s_val==3
        F=delta_r* ((C^alpha)^(1-psi)) /(1-psi);
    end
end

if s_val==4 % If HH is dead, just overrule everything and give zero utility
    F=0;
end

end