% Loads US male data needed by HuggettVentura1999:
%  HV1999_sj    conditional survival prob, ages 20..100  (length 81)
%  HV1999_ybar  age-earnings profile e_bar_j, ages 20..100 (length 81)
%
% Survival: HV1999 Section 3.1 cites the SSA 1994 male life table (footnote
% 15: thanks to J. Gokhale). We use Bell, Wade & Goss (1992, Actuarial
% Study 107) Table 7 male death probs (q_x at exact ages 0,30,60,65,70,100
% per decade 1900..2080), interpolate to all ages 0..100 and to year 1994.
%
% Age-earnings: Hansen (1993) Table II male weights at age-group centres,
% linearly interpolated and taken to zero at age 75, exactly as HV1999 fn16
% describes. See the block below.

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

% Interpolate to all ages 0..100 (rows = years, cols = ages 0..100).
% LOG-linear in age, not linear. Table 7 gives q_x only at 0, 30, 60, 65, 70 and
% 100, so the 0->30 and 30->60 gaps are being bridged across 30 years. Adult
% mortality is approximately Gompertz, q_x ~ exp(a+b*x), so q_x is convex in age
% and linear interpolation OVERSTATES it between anchors -- which made survival
% too low: cumulative survival 20->65 came out 0.687 against roughly 0.76 for
% U.S. males in 1994, and the model's 65-100 population share was 0.1345 against
% roughly 0.17-0.18 in the data. Too few agents reaching retirement means too
% little total social security spending, hence too low an equilibrium theta,
% hence too little labour distortion and too small a system for the Boskin
% reform to act on -- the two residual gaps as of 2026-10-06.
% Interpolating log(q_x) and exponentiating gives 20->65 survival of 0.732 and a
% retiree share of 0.1591, and 20->80 survival of 0.372 against about 0.37 in
% the data. Still short of the true per-age SSA 1994 male life table, which is
% what should eventually replace this.
BWG_deathprobs_age=exp(interp1(BWG_ExactAge',log(BWG_Male'),(0:1:100)')');
% Interpolate to all years 1900..2080
BWG_deathprobs=interp1(BWG_Year,BWG_deathprobs_age,(1900:1:2080));

% HV1999 uses US 1994 -> row index 95. Ages 20..100 -> column indexes 21..101.
HV1999_sj=1-BWG_deathprobs(95,21:101);
HV1999_sj=HV1999_sj(:); % column
HV1999_sj(end)=0; % death certain at terminal age 100

clear BWG_Year BWG_ExactAge BWG_Male BWG_deathprobs_age BWG_deathprobs

%% Age-earnings profile: Hansen (1993) raw weights, built as HV1999 fn16 states
% HV1999 footnote 16: "Hansen estimates median wage rates in cross-sectional
% data for males in different age groups. We use his values for the center of
% the age group and linearly interpolate to get the remaining values. We set the
% labor endowment in efficiency units of agents at a real-life age of 75 to
% zero." Fig 1 (p.511) is the resulting profile.
%
% Hansen (1993) Table II, "Weights assigned to age-sex", MALE column (the female
% column is 0.52,0.69,0.89,0.90,0.87,0.84,0.66 and is not used here). Sourced
% from the ConesaKrueger1999 replication in vfitoolkit-matlab-replication, which
% quotes the table directly:
%     16-19  0.56     45-54  1.39
%     20-24  0.78     55-64  1.33
%     25-34  1.14     65+    0.89
%     35-44  1.37
% Placed at age-group CENTRES and interpolated linearly, per fn16.
%
% The 65+ weight is NOT used: it is an open-ended group with no centre, and
% interpolating 59.5 -> a 65+ centre would give a gentle decline, whereas Fig 1
% plunges after 60 and reaches zero at 75. Running straight from the 55-64
% centre down to zero at 75 reproduces Fig 1: it gives 0.944 at age 64 against
% 0.940 read off the figure. It is also what fn16's "zero at 75" instruction
% implies. Ages 65+ never work anyway (Jr=46 is real age 65), so only the
% 20-64 part of this profile enters the model.
%
% This REPLACES a visual digitisation of Fig 1 (2026-10-05..06), which had the
% peak at age 44 rather than 49 and ran 9-14% low through the rise and plateau.
% The level difference is largely harmless -- a uniform scaling of kappa_j is
% absorbed by w in general equilibrium, and now also by ybar_mean_unwage, which
% is itself a GE unknown -- but the peak location is a real shape error.
HV1999_Hansen1993_centre=[17.5, 22, 29.5, 39.5, 49.5, 59.5];
HV1999_Hansen1993_male  =[0.56, 0.78, 1.14, 1.37, 1.39, 1.33];
HV1999_ybar_ages=(20:100)';
HV1999_ybar=interp1(HV1999_Hansen1993_centre,HV1999_Hansen1993_male,HV1999_ybar_ages,'linear','extrap');
% From the last centre (55-64) run linearly down to zero at 75, and stay zero.
HV1999_ybar_tail=HV1999_ybar_ages>=HV1999_Hansen1993_centre(end);
HV1999_ybar(HV1999_ybar_tail)=HV1999_Hansen1993_male(end) ...
    *max(0,(75-HV1999_ybar_ages(HV1999_ybar_tail))/(75-HV1999_Hansen1993_centre(end)));
HV1999_ybar(HV1999_ybar_ages>=75)=0;
% HV1999_Hansen1993_centre/_male are deliberately NOT cleared: they are Hansen's
% own tabulated values and Fig 1 plots them as markers on top of the
% interpolated profile. Clearing them made that figure silently drop its data
% source -- the plotting code guards with exist(), so it just drew the line.
% Six doubles each; keeping them costs nothing.
clear HV1999_ybar_ages HV1999_ybar_tail
