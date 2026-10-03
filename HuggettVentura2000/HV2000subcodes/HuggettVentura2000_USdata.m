% Loads US 1990 male data needed by HuggettVentura2000:
%  HV2000_sj    conditional survival prob, ages 20..100  (length 81)
%  HV2000_ybar  age-earnings profile e_bar_j, ages 20..100 (length 81)
%
% Survival: HV2000 Section 4.1 cites Social Security Administration (1992),
% which is Actuarial Study 107 (Bell, Wade & Goss 1992). Table 7 of BWG1992
% reports period probabilities of death within one year (q_x) at exact ages
% 0, 30, 60, 65, 70, 100 for each decade 1900..2080. We interpolate to
% all ages 0..100 and to year 1990.
%
% Age-earnings: HV2000 Fig 2 is built from US 1990 male median earnings
% (Social Security Bulletin 1995 Table 4.B6) times labor-force participation
% rates (Fullerton 1992). Pending exact-source recovery, this script uses
% the Conesa-Krueger (1999) / DIS1999 working-age profile, rescaled so the
% peak matches HV2000 Fig 2 (~1.55), then zero from retirement age onward.

%% Bell, Wade & Goss (1992) Table 7 -- male death probs at selected exact ages
% Rows: year 1900,1910,...,2080. Cols: ages 0, 30, 60, 65, 70, 100.
BWG_Year=(1900:10:2080)';
BWG_ExactAge=[0,30,60,65,70,100];
BWG_Male=[0.14595, 0.00838, 0.0293,  0.04158, 0.06182, 0.44999;...
          0.12006, 0.00713, 0.02968, 0.04151, 0.06173, 0.45274;...
          0.08593, 0.0066,  0.02481, 0.03714, 0.05728, 0.43876;...
          0.06495, 0.00491, 0.0274,  0.03945, 0.05751, 0.4475;...
          0.05286, 0.0034,  0.02663, 0.03812, 0.05611, 0.44203;...
          0.03279, 0.00213, 0.02476, 0.03487, 0.05046, 0.38648;...
          0.02937, 0.00183, 0.02392, 0.03515, 0.05019, 0.38224;...
          0.02246, 0.00209, 0.02348, 0.03416, 0.04887, 0.36256;...
          0.01398, 0.00189, 0.01843, 0.02881, 0.04312, 0.34225;...
          0.01019, 0.00212, 0.01562, 0.02479, 0.03763, 0.34199;...
          0.00723, 0.00238, 0.01335, 0.02191, 0.03418, 0.3349;...
          0.00603, 0.00184, 0.01189, 0.02009, 0.03195, 0.32487;...
          0.00551, 0.00178, 0.01115, 0.01898, 0.03025, 0.3075;...
          0.00507, 0.00174, 0.0105,  0.01798, 0.0287,  0.29138;...
          0.00469, 0.0017,  0.0099,  0.01708, 0.02729, 0.27655;...
          0.00435, 0.00166, 0.00936, 0.01625, 0.02601, 0.26294;...
          0.00405, 0.00162, 0.00887, 0.0155,  0.02484, 0.25042;...
          0.00378, 0.00159, 0.00842, 0.01481, 0.02377, 0.2389;...
          0.00355, 0.00155, 0.008,   0.01417, 0.02278, 0.22828];

% Interpolate to all ages 0..100 (rows = years, cols = ages 0..100)
BWG_deathprobs_age=interp1(BWG_ExactAge',BWG_Male',(0:1:100)')';
% Interpolate to all years 1900..2080
BWG_deathprobs=interp1(BWG_Year,BWG_deathprobs_age,(1900:1:2080));

% HV2000 uses US 1990 -> row index 91. Ages 20..100 -> column indexes 21..101.
HV2000_sj=1-BWG_deathprobs(91,21:101);
HV2000_sj=HV2000_sj(:); % column
HV2000_sj(end)=0; % death certain at terminal age 100

clear BWG_Year BWG_ExactAge BWG_Male BWG_deathprobs_age BWG_deathprobs

%% Age-earnings profile (HV2000 Figure 2)
% ybar_j = (1990 male median earnings of workers) × (1990 male LFP rate), by age bin.
% HV2000 cite the SSA Annual Statistical Supplement (1995, Table 4.B6) for the median
% earnings and Fullerton 1992 MLR for the LFP rates; both underlying series are public.

% Midpoint (real-life age) of each SSA Table 4.B6 bin. Used for linear interpolation
% to single-year ages 20-64.
SSA_bin_mid = [...
    18,   ... % "Under 20" (SSA reports the 16-19 group; midpt 17.5 rounded to 18)
    22,   ... % "20-24"
    27,   ... % "25-29"
    32,   ... % "30-34"
    37,   ... % "35-39"
    42,   ... % "40-44"
    47,   ... % "45-49"
    52,   ... % "50-54"
    57,   ... % "55-59"
    60.5, ... % "60-61"
    63,   ... % "62-64"
    67,   ... % "65-69"
    70.5, ... % "70-71"
    77    ... % "72 or older" (representative age picked; bins past 64 don't enter the model)
    ];

% Median earnings of MALE workers in 1990 (current US dollars).
% Source: SSA Annual Statistical Supplement 2003, Table 4.B6, "Men" panel, 1990 row.
% Downloaded from https://www.ssa.gov/policy/docs/statcomps/supplement/2003/4b.pdf
% (saved as 4b.pdf in this folder). Verified row-by-row against the PDF.
SSA_med_earn = [...
    2058,  ... % Under 20
    8945,  ... % 20-24
    16412, ... % 25-29
    21211, ... % 30-34
    24424, ... % 35-39
    27608, ... % 40-44
    29074, ... % 45-49 (peak)
    28207, ... % 50-54
    25509, ... % 55-59
    23243, ... % 60-61
    17408, ... % 62-64
    7714,  ... % 65-69
    6153,  ... % 70-71
    5129   ... % 72+
    ];

% Male civilian labor-force participation rates, 1990 US annual averages,
% by SSA age bin. Source: Fullerton, H.N. Jr. (1992) "Evaluation of labor
% force projections to 1990", Monthly Labor Review August 1992, Table 1,
% "Actual 1990" column, Men panel. Same article HV2000 cite (HV2000 fn 13).
% PDF: https://www.bls.gov/opub/mlr/1992/08/art1full.pdf
% (also saved as Fullerton1992.pdf in this folder).
% Fullerton Table 1 reports 25-34, 35-44, 45-54, 60-64, and 70+ as aggregates
% (no 5-year subbins); the aggregate value is reused for each SSA subbin that
% falls within. The SSA 60-61 vs 62-64 split (early-SS-retirement effect) and
% the 70-71 vs 72+ split are not published by Fullerton and remain flagged
% TEMPORARY -- approximated to be consistent with Fullerton's 60-64 = 0.555
% (pop-weighted 0.685, 0.465 at 2:3 = 0.553) and 70+ = 0.108.
BLS_LFP_approx = [...
    0.557, ... % Under 20 (16-19; pop-weighted avg of Fullerton 16-17=0.437 and 18-19=0.670)
    0.843, ... % 20-24 (Fullerton, exact)
    0.942, ... % 25-29 (Fullerton 25-34 aggregate; subbin not published)
    0.942, ... % 30-34 (Fullerton 25-34 aggregate; subbin not published)
    0.944, ... % 35-39 (Fullerton 35-44 aggregate; subbin not published)
    0.944, ... % 40-44 (Fullerton 35-44 aggregate; subbin not published)
    0.907, ... % 45-49 (Fullerton 45-54 aggregate; subbin not published)
    0.907, ... % 50-54 (Fullerton 45-54 aggregate; subbin not published)
    0.798, ... % 55-59 (Fullerton, exact)
    0.685, ... % 60-61 - TEMPORARY; Fullerton aggregates 60-64=0.555
    0.465, ... % 62-64 - TEMPORARY; same caveat (sharp early-SS-eligibility drop)
    0.260, ... % 65-69 (Fullerton, exact)
    0.108, ... % 70-71 (Fullerton 70+ aggregate)
    0.108  ... % 72+   (Fullerton 70+ aggregate)
    ];

realage_working=20:64;
ybar_working=interp1(SSA_bin_mid,SSA_med_earn.*BLS_LFP_approx,realage_working,'linear','extrap');
ybar_working=ybar_working/max(ybar_working)*1.55; % rescale so peak == 1.55, matching HV2000 Fig 2 visually

HV2000_ybar=zeros(81,1);
HV2000_ybar(1:45)=ybar_working;
clear SSA_bin_mid SSA_med_earn BLS_LFP_approx realage_working ybar_working
