function F=Aiyagari1994_ReturnFn(aprime_val, a_val, s_val,alpha,delta,mu,r)
% You can include comments in the return function file, but only on
% seperate lines to the commands. If you try putting a comment on the same
% line as a command (and following line is not blank) then this will cause
% an error with the GPU codes.

F=-Inf;
w=(1-alpha)*((r+delta)/alpha)^(alpha/(alpha-1));
c=w*s_val+(1+r)*a_val-aprime_val; 
%c=wl_t+(1+r)a_t-a_{t+1}
if c>0
    if mu==1
        F=log(c);
    else
        F=(c^(1-mu) -1)/(1-mu);
    end
end

end