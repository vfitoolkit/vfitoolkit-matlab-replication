function sr=HV2000_SavingRateFn(d,aprime,a,ebar,z,agej,Jr,r,w,T,b_common,ybar_mean,ybar_j,g,haveSS)
% Per-agent saving rate = Savings / Income. Used by Figure 7 to take the age-conditional MEDIAN
% (the paper's preferred summary; AllStats.SavingRate.Mean would be the population avg, not the same thing).
% Savings = (1+g)*aprime - a (growth-adjusted change in assets). Income mirrors HV2000_IncomeFn.

savings=(1+g)*aprime-a;

if agej<Jr
    earnings=z*ybar_j*w;
    benefit=0;
else
    earnings=0;
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

income=earnings+r*a+T+benefit;

if income>0
    sr=savings/income;
else
    sr=0;
end

end
