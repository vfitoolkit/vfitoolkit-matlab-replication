function hist_idx_prime=LaunWallenius2015_aprimeFn(d,hist_idx,agej,agej_pb_first,agej_di_max)
% Experienceasset law of motion for the history-index state.
%
% State coordinates encoded in hist_idx (see driver header for layout):
%   swa = first stop-work age (model period; 0 = never stopped)
%   dia = first DI-claim age  (model period; 0 = never claimed DI)
%   pba = first PB-claim age  (model period; 0 = never claimed PB)
%
% Action codes d:
%   1: keep working, no claim change
%   2: stop work (or stay stopped), no claim change
%   3: apply for DI (forces stop work)
%   4: claim PB and stop working
%   5: claim PB and keep working
% Eligibility constraints (h<=2 for DI, etc.) are enforced in ReturnFn via -Inf.
% This aprimeFn just computes the formal transition; infeasible (d, state)
% combinations get their utility set to -Inf in the ReturnFn so the policy
% function never lands on them anyway.

NEVER=0;
n_pba_options=21;          % {NEVER, 37..56}
n_no_di=57*n_pba_options;  % 1197 indices for the no-DI branch

% --- Decode hist_idx -> (swa, dia, pba) ---
if hist_idx<=n_no_di
    swa=floor((hist_idx-1)/n_pba_options);              % 0..56
    pba_idx=mod(hist_idx-1,n_pba_options)+1;            % 1..21
    if pba_idx==1
        pba=NEVER;
    else
        pba=pba_idx+35;                                 % 37..56
    end
    dia=NEVER;
else
    dia=hist_idx-n_no_di;                               % 1..40
    swa=dia;                                            % work absorbing on DI
    pba=NEVER;
end

% --- Apply d-dispatch (with eligibility no-ops when infeasible) ---
swa_new=swa;
dia_new=dia;
pba_new=pba;

if d==1
    % keep working, no claim change; only feasible if swa==NEVER (ReturnFn enforces)
elseif d==2
    if swa==NEVER
        swa_new=agej;
    end
elseif d==3
    if dia==NEVER && agej<=agej_di_max && pba==NEVER
        dia_new=agej;
        swa_new=agej;
    end
elseif d==4
    if pba==NEVER && agej>=agej_pb_first && dia==NEVER
        pba_new=agej;
        if swa==NEVER
            swa_new=agej;
        end
    end
elseif d==5
    if pba==NEVER && agej>=agej_pb_first && swa==NEVER && dia==NEVER
        pba_new=agej;
    end
end

% --- Encode (swa_new, dia_new, pba_new) -> hist_idx_prime ---
if dia_new>NEVER
    hist_idx_prime=n_no_di+dia_new;
else
    if pba_new==NEVER
        pba_idx_new=1;
    else
        pba_idx_new=pba_new-35;                          % 2..21
    end
    hist_idx_prime=swa_new*n_pba_options+pba_idx_new;
end

end
