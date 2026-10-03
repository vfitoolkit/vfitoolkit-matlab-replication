function F=HuggettVentura2000_ReturnFn(d,aprime,a,ebar,z,agej,Jr,sigma,g,r,w,tau,theta,ybar_j,T,b_common,ybar_mean,borrowFlag,haveSS)
% HV2000 budget constraint (eqn 5'), CRRA period utility (eqn 2).
% d is singleton dummy; z is productivity shock (=1 in Model 1).
% Bend points and AIME cap are derived from ybar_mean (Section 4.3):
% b_bend1=0.20*ybar_mean, b_bend2=1.24*ybar_mean, bend rates 90/32/15%.

F=-Inf;

% Labor earnings (zero after retirement)
if agej<Jr
    earnings=z*ybar_j*w;
else
    earnings=0;
end

% Social security benefit (eqn 7)
if agej<Jr
    benefit=0;
else
    b_bend1=0.20*w*ybar_mean;
    b_bend2=1.24*w*ybar_mean;
    if ebar<=b_bend1
        b_er=0.90*ebar;
    elseif ebar<=b_bend2
        b_er=0.90*b_bend1+0.32*(ebar-b_bend1);
    else
        b_er=0.90*b_bend1+0.32*(b_bend2-b_bend1)+0.15*(ebar-b_bend2);
    end
    benefit=haveSS*(b_common+b_er/((1+g)^(agej-Jr)));
end

% Budget: c + (1+g)*aprime = a*(1+r*(1-tau)) + (1-tau-theta)*earnings + T + benefit
c=a*(1+r*(1-tau))+(1-tau-theta)*earnings+T+benefit-(1+g)*aprime;

if c>0
    F=(c^(1-sigma))/(1-sigma);
end

if aprime<-borrowFlag*w % a̲ = -borrowFlag*w; 0 or -ŵ per Params.borrowFlag (HV2000 Table 2)
    F=-Inf;
end

end
