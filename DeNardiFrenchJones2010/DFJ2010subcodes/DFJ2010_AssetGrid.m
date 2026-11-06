% DFJ2010_AssetGrid: the asset grid, shared by DeNardiFrenchJones2010.m and DFJ2010_DataWork.m (the initial
% distribution is placed on this grid, so the two must agree; the model script checks this when it loads the data).
%
% Piecewise grid (2026-09-23): the GMM targets are median assets by cohort x PI quintile x year, and every one of
% the 102 target medians lies between 0 and 280,000 (1998 dollars), while the previous cubic grid on [0, 2m]
% (n_a=501) put 235 of its points above 300,000. Medians sit on grid points, so the objective is flat in a
% parameter until the median crosses a grid point; the coarse steps at 100-300k (1.6-3.4k) were why the bequest
% parameters theta and k did not move in the estimation. Now 801 points on [0, 300k] (quadratic, denser near zero
% where the value fn has most curvature: steps of 0.14k at 10k, 0.43k at 100k, 0.75k at 300k) and 200 points on
% (300k, 2m] (power 1.5, steps from 0.6k to 13k; only the value fn uses this range, the targets never do).
% The upper bound 2m is DFJ2010's asset top-code (asstmax in getdat25_wave8.g).
n_a=1001;
a_grid=[300000*linspace(0,1,801)'.^2; 300000+1700000*linspace(0,1,201)'.^1.5];
a_grid=a_grid([1:801,803:end]); % the lower segment ends at 300k and the upper begins there; drop the duplicate
