function F=Huggett1993_ReturnFn(aprime_val, a_val, z_val,mu,q)
% You can include comments in the return function file, but only on
% seperate lines to the commands. If you try putting a comment on the same
% line as a command (and following line is not blank) then this will cause
% an error with the GPU codes.

% Huggett (1993) refers to mu as sigma.

F=-Inf;
c=a_val+z_val-q*aprime_val; 
%c=a_t+e_t-q*a_{t+1}
if c>0
    if mu==1
        F=log(c);
    else
        F=(c^(1-mu) -1)/(1-mu);
    end
end

end