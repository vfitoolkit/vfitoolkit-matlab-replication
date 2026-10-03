function ebarprime=HuggettVentura2000_ebarprimeFn(d,ebar,z,agej,Jr,w,ybar_j,ybar_mean,g)
% Law of motion for average past earnings (HV2000 eqn 8).
% d is a singleton dummy (HV2000 has no labor choice; required by experienceassetz signature).
% z is the productivity shock (singleton=1 in Model 1; meaningful in Models 2-4).
% AIME-input cap ebar_max = 2.47*ybar_mean (Section 4.3); ybar_mean is the
% population-average earnings, endogenous in GE.

ebarprime=0; % keep GPU happy

if agej<Jr
    ebar_max=2.47*w*ybar_mean;
    earnings=min(z*ybar_j*w,ebar_max);
    ebarprime=(ebar*(agej-1)+earnings)/agej;
else
    ebarprime=ebar/(1+g);
end

end
