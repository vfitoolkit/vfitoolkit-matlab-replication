function ebarprime=HuggettVentura1999_ebarprimeFn(l,ebar,logprod,agej,Jr,kappa_j,ybar_mean_unwage)
% logprod is LOG productivity, and it occupies the slot that is a per-type
% PARAMETER in the no-idiosyncratic spec (logabilityi, ability fixed at birth)
% and an exogenous STATE in the idiosyncratic spec (z, the AR(1) deviation from
% the age-mean log endowment). Same position, same meaning, so one file serves
% both; only the driver's anonymous-function wrapper differs.
% HV1999 eqn (7), rewritten in "ratio to mean working-age earnings" units:
%   ebar (state) := AIME / (w * ybar_mean_unwage)
% so the grid points {0, 0.20, 1.24, 2.47} match the bend points / AIME cap
% directly and never need rebuilding when the GE wage moves.
%
% Dividing eqn (7) by (w*ybar_mean_unwage) and noting that w cancels:
%   ebar' = [ebar*(agej-1) + min(e(z,j)*l/ybar_mean_unwage, 2.47)] / agej   (j<Jr)
%   ebar' = ebar                                                            (j>=Jr)
% After retirement the (1+g)^{j-R} detrend stays inside the benefit formula
% in the ReturnFn (HV1999 eqn 15).

ebarprime=ebar;

if agej<Jr
    earnings_ratio=min(exp(logprod)*kappa_j*l/ybar_mean_unwage,2.47);
    ebarprime=(ebar*(agej-1)+earnings_ratio)/agej;
end

end
