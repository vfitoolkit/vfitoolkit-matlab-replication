% Guvenen2007.m  --  MASTER DRIVER for the Guvenen (2007 AER) baseline replication.
% Set doPart to choose which stages to run; each stage writes its own diary/results/figures
% and can also be run on its own. Data work first (CPU), then the model pipeline (GPU), then
% the delta/grid diagnostics.
%
%   doPart(1): data work        Guvenen2007_DataWork   -> the empirical var(log c) age
%              profile + 95% bootstrap CIs (CEX Interview, Krueger-Perri nondurables,
%              Deaton-Paxson age-cohort decomposition); writes SavedOutput/Guvenen2007_DataWork.mat
%   doPart(2): RunRIP           RIP model (plumbing check)          RunRIP_results.mat
%   doPart(3): RunBaseline      HIP model, lambda in {0,.4,.62,1}   RunBaseline_results.mat
%   doPart(4): RunEducationGroups  college/HS x HIP/RIP (Fig 11)    RunEducationGroups_results.mat
%   doPart(5): MakeReplicationTable  table + Fig 7/8 + Fig 11       replication_table.tex, figs
%   doPart(6): RunBaseline_deltadiag delta=0.966 fixed             _deltadiag diary/results
%   doPart(7): RunBaseline_gridtest  wealth-grid sweep @0.966      _gridtest diary/results
%   doPart(8): RunBaseline_gridconv  rise convergence over n_a     _gridconv diary/results
%
% NOTE doPart(1) is CPU-only and can be run in any MATLAB; doPart(2)-(4),(6)-(8) solve the
% model and need the GPU. doPart(5) is CPU post-processing of the saved .mat results.
% Each sub-script starts with `clearvars -except doPart` so doPart survives across stages.
%
% The stage scripts and their helpers (BuildPiZ, KalmanSetup, ReturnFn_HIP/RIP) live in
% ./G2007_subcodes/; MATLAB does not search subfolders of the working directory, so we
% addpath it here. Run this master from baseline/ so the './SavedOutput/...', '../../DataSets/...'
% etc. relative paths inside the subcodes still resolve to here.
addpath('./G2007_subcodes');

doPart=[0,1,1,1,1,0,0,0]; % (1) data work cached (SavedOutput/Guvenen2007_DataWork.mat); (6)-(8) delta/grid diagnostics already done (see NOTES.md / tex Appendix B) -- set any back to 1 only to rebuild/rerun

fprintf('\n################ Guvenen (2007) replication -- full pipeline ################\n');
fprintf('start %s   doPart=[%s]\n',char(datetime('now')),num2str(doPart));

if doPart(1)==1
    fprintf('\n>>> [1/8] Guvenen2007_DataWork  (%s)\n',char(datetime('now')));
    Guvenen2007_DataWork
end
if doPart(2)==1
    fprintf('\n>>> [2/8] RunRIP  (%s)\n',char(datetime('now')));
    RunRIP
end
if doPart(3)==1
    fprintf('\n>>> [3/8] RunBaseline  (%s)\n',char(datetime('now')));
    RunBaseline
end
if doPart(4)==1
    fprintf('\n>>> [4/8] RunEducationGroups  (%s)\n',char(datetime('now')));
    RunEducationGroups
end
if doPart(5)==1
    fprintf('\n>>> [5/8] MakeReplicationTable  (%s)\n',char(datetime('now')));
    MakeReplicationTable
end
if doPart(6)==1
    fprintf('\n>>> [6/8] RunBaseline_deltadiag  (%s)\n',char(datetime('now')));
    RunBaseline_deltadiag
end
if doPart(7)==1
    fprintf('\n>>> [7/8] RunBaseline_gridtest  (%s)\n',char(datetime('now')));
    RunBaseline_gridtest
end
if doPart(8)==1
    fprintf('\n>>> [8/8] RunBaseline_gridconv  (%s)\n',char(datetime('now')));
    RunBaseline_gridconv
end

fprintf('\n################ Guvenen (2007) pipeline finished %s ################\n',char(datetime('now')));
