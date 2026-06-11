function ap=LaunWallenius2015_AP(swa,dia,agej,agejshifter,w_c0_i,w_c1_i,w_c2_i,lbar,ppmax,N_j)
% Average Pension Points: top-15-best rule with 3-year-pre-DI projection
% (LW2015 p.130-131).
%
% Pure scalar implementation matching the OPB scalar-loop style: no local
% arrays, no sort. Closed-form wage and pension-points computed inline.
% Top-15 found via 15 passes through the per-age pp values, with each pass
% selecting the next-largest (pp, j) in lexicographic order (pp descending,
% j ascending on ties). The j-tiebreak handles exact ties (e.g. ages capped
% at ppmax, or zero-pp ages).
%
% Cost: 15 * N_j scalar pp evaluations + lex comparisons. ~25x more work
% than a sort-based vectorized version, but stays within the scalar-only
% GPU-arrayfun envelope used by the OPB scalar refactor.

NEVER=0;
agej_di_max=40;  % age 64 -> model age 40

% --- Projection pp (computed once; scalar accumulator over up-to-3 ages) ---
if dia>NEVER
    proj_sum=0;
    n_proj=0;
    for j_proj=max(1,dia-3):dia-1
        real_age_pj=j_proj+agejshifter;
        w_pj=w_c2_i*real_age_pj^2+w_c1_i*real_age_pj+w_c0_i;
        proj_sum=proj_sum+w_pj;
        n_proj=n_proj+1;
    end
    if n_proj>0
        w_proj=proj_sum/n_proj;
    else
        w_proj=0;
    end
    pp_proj=max(0,min(w_proj*lbar,ppmax+1)-1);
else
    pp_proj=0;
end

% --- Boundary of "real working" ages ---
if dia>NEVER
    last_real_working_j=dia-1;
elseif swa>NEVER
    last_real_working_j=swa-1;
else
    last_real_working_j=agej-1;
end

% --- Top-15 via 15 passes, lexicographic next-largest selection ---
top_sum=0;
prev_pp=Inf;
prev_j=0;
for pass_k=1:15
    best_pp=-1;
    best_j=0;
    for j=1:N_j
        % pp at age j (closed-form, scalar)
        if j<=last_real_working_j
            real_age_j=j+agejshifter;
            w_j=w_c2_i*real_age_j^2+w_c1_i*real_age_j+w_c0_i;
            pp_j=max(0,min(w_j*lbar,ppmax+1)-1);
        elseif (dia>NEVER) && (j>=dia) && (j<=agej_di_max)
            pp_j=pp_proj;
        else
            pp_j=0;
        end

        % Eligible: strictly less than prev (pp, j) in lex order
        % i.e. pp_j < prev_pp, OR (pp_j == prev_pp AND j > prev_j)
        eligible=(pp_j<prev_pp)||((pp_j==prev_pp)&&(j>prev_j));
        if eligible
            % Beats current best in lex order (pp desc, j asc on ties)?
            if (pp_j>best_pp)||((pp_j==best_pp)&&(j<best_j))
                best_pp=pp_j;
                best_j=j;
            end
        end
    end
    top_sum=top_sum+best_pp;
    prev_pp=best_pp;
    prev_j=best_j;
end

ap=top_sum/15;

end
