function F=HopenhaynRogerson1993_ReturnFn(aprime_val, a_val, z_val, p, alpha, tau, cf)
% a_val is what HopenhaynRogerson1993 call n_{t-1}, aprime_val is n_t.

F=-Inf;

gnnlag=0; % pg 918 g(n_t,n_{t-1}) in notation of Hopenhayn & Rogerson (1993), the labour adjustment-costs (firing cost)
if aprime_val<a_val
    gnnlag=-tau*(aprime_val-a_val); % Note that gnnlag>=0, it is 'minus a negative number'
end

% pg 918
pi=p*z_val*(aprime_val^alpha)-aprime_val-p*cf-gnnlag; %Hardcodes normalization w=1

% Following few lines are just required to deal with Footnote 5.
if a_val==10^6 % Impose Footnote 5 (see first few lines of HopenhaynRogerson1993 for explanation of this)
    gnnlag=0; % The 'actual' number of lag employees of a new entrant is zero, so firing costs must be zero.
    pi=p*z_val*(aprime_val^alpha)-aprime_val-gnnlag; % same as above, but without the -p*cf term.
end
if aprime_val==10^6 % Impose Footnote 5 (see first few lines of HopenhaynRogerson1993 for explanation of this)
    pi=-Inf; % Existing firms cannot choose to become new entrants.
end

F=pi;

end