function F=HV2000_ReturnFn_M4(d,aprime,a,ebar,z,e,agej,Jr,sigma,g,r,w,tau,theta,ybar_j,T,b_common,ybar_mean,borrowFlag,haveSS)
% Model 4 ReturnFn: same as Models 1-3 plus the iid e shock multiplying
% earnings. e is a level multiplier exp(z_2) (HV2000 footnote 12).

F=-Inf;

if agej<Jr
    earnings=z*e*ybar_j*w;
else
    earnings=0;
end

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

c=a*(1+r*(1-tau))+(1-tau-theta)*earnings+T+benefit-(1+g)*aprime;

if c>0
    F=(c^(1-sigma))/(1-sigma);
end

if aprime<-borrowFlag*w % a̲ = -borrowFlag*w; 0 or -ŵ per Params.borrowFlag (HV2000 Table 2)
    F=-Inf;
end

end
