function F=Hansen1985_ReturnFn(d_val, aprime_val, a_val, z_val, alpha, delta, A, B, Economy)

F=-Inf;
C=z_val*(a_val^alpha)*(d_val^(1-alpha))-(aprime_val-(1-delta)*a_val);
if C>0
    if Economy==1
        F=log(C)+A*log(1-d_val); % Divisible labour
    elseif Economy==2
        F=log(C)+B*(1-d_val); % Indivisible labour
    end
end

end