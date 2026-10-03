function ssb=HV2000_SSBenefitFn(d,aprime,a,ebar,z,agej,Jr,g,b_common,ybar_mean,w,haveSS)
% Social security benefit per agent (HV2000 eqn 7). Used by FnsToEvaluate
% as the per-agent benefit so that AggVars gives total SS outlays per
% population (then multiplied by appropriate weights via mewj).

ssb=0;

if agej>=Jr
    b_bend1=0.20*w*ybar_mean;
    b_bend2=1.24*w*ybar_mean;
    if ebar<=b_bend1
        b_er=0.90*ebar;
    elseif ebar<=b_bend2
        b_er=0.90*b_bend1+0.32*(ebar-b_bend1);
    else
        b_er=0.90*b_bend1+0.32*(b_bend2-b_bend1)+0.15*(ebar-b_bend2);
    end
    ssb=haveSS*(b_common+b_er/((1+g)^(agej-Jr)));
end

end
