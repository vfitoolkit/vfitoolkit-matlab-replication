function TWpct = HV2000_KotlikoffSummersStat(Params, K)
% Kotlikoff-Summers (1981) transfer-wealth-as-share-of-total-wealth statistic.
%
% HV2000 Section 5.1 / Table 4 report the "Transfer wealth (%)" column using this
% decomposition. Every alive agent in HV2000 receives a lump-sum transfer T each
% period; the current-value (compounded at after-tax rate) of all past T-receipts
% is that agent's transfer wealth.
%
% Growth-adjusted per-period return: Rhat = (1 + r*(1-tau)) / (1 + g)
% Per-age transfer wealth: TW_j = T * sum_{k=0}^{j-1} Rhat^k = T*(Rhat^j-1)/(Rhat-1)
%                         (=  T*j if Rhat == 1)
% Aggregate TW = sum_j mewj * TW_j; percentage of total wealth = 100 * TW / K.
%
% For the PType variant (Params.T is a vector of length N_i), aggregate TW is
% linear in T so we can just replace T with the ptype-mass-weighted mean.
%
% Inputs
%   Params  Struct with fields r, tau, g, T, mewj, J
%           (if PType: also alphai_dist matching length(T))
%   K       Aggregate capital (typically AllStats.K.Mean)
%
% Output
%   TWpct   Transfer wealth as a percentage of aggregate capital

if isscalar(Params.T)
    T_bar=Params.T;
else
    % PType case: aggregate T weighted by ptype distribution
    T_bar=sum(Params.alphai_dist(:).*Params.T(:));
end

Rhat=(1+Params.r*(1-Params.tau))/(1+Params.g);
if abs(Rhat-1)>1e-10
    TW_j=T_bar*(Rhat.^(1:Params.J)-1)/(Rhat-1);
else
    TW_j=T_bar*(1:Params.J);
end

TWpct=100*sum(Params.mewj(:)'.*TW_j)/K;

end
