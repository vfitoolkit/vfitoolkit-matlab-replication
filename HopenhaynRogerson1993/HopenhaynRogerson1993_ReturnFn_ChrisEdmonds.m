function F=HopenhaynRogerson1993_ReturnFn_ChrisEdmonds(aprime_val, a_val, z_val, p, alpha, tau, cf)
% a_val is what HopenhaynRogerson1993 call n_{t-1}, aprime_val is n_t.

F=-Inf;

gnnlag=0; % pg 918 g(n_t,n_{t-1}) in notation of Hopenhayn & Rogerson (1993), the labour adjustment-costs (firing cost)
if aprime_val<a_val
    gnnlag=-tau*(aprime_val-a_val); % Note that gnnlag>=0, it is 'minus a negative number'
end

% pg 918
pi=p*z_val*(aprime_val^alpha)-aprime_val-cf-gnnlag; % Only difference is -cf, instead of -p*cf

F=pi;

end