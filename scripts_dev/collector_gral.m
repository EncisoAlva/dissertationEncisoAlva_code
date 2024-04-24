% This script is paired with
%       generator04
% It is meant to search for all test cases, solve each one with all
% available solver, and then run all evaluation metrics. The results are
% first saved on one file per test case (for robustness) and then condensed
% on a big file to compare the methods.
%
% Author: Julio C Enciso-Alva (2023)
%         juliocesar.encisoalva@mavs.uta.edu
%

cd(fileparts(which(mfilename)));
digits(60)

% which test data will be loaded
BaseName   = 'protocol_04_test2';
hard_reset = false;

%% SETUP
originalPath = pwd;
addpath(pwd)
cd('..')
basePath  = pwd;

%% INFO ABOUT TEST CASES
cd(['./data/', BaseName])
dataPath = pwd;

load("metadata.mat")
load("metadata2.mat")

nCases      = info.nTrials;
nConditions = length(info.SNRvals);
ConditionFolder = cell(nConditions,1);
for i = 1:nConditions
  ConditionFolder{i} = [ 'SNR_', num2str(info.SNRvals(i)) ];
end

cd(basePath)

%% LOAD SOLVERS
cd('solvers')
lst = dir('*.m');

nSolvers   = length(lst);
SolveNames = cell(nSolvers,1);
for i = 1:nSolvers
  tmp = split(lst(i).name,'.');
  SolveNames{i} = tmp{1};
end

cd('..');

%% LOAD EVALUATION METRICS
cd('metrics')
lst = dir('*.m');

nEvals    = length(lst);
EvalNames = cell(nEvals+3,1);
for i = 1:nEvals
  tmp = split(lst(i).name,'.');
  EvalNames{i} = tmp{1};
end
% additional variables:
EvalNames{nEvals+1} = 'Run_Time';
%EvalNames{nEvals+2} = 'ParTune_Time';
EvalNames{nEvals+2} = 'Depth';
EvalNames{nEvals+3} = 'Kappa';

cd('..');

% container variable
tableNames = cell(1, nEvals+3+3 );
for ii = 1:(nEvals+3)
  tableNames{ii} = EvalNames{ii};
end
tableNames{nEvals+3 +1}  = 'SNR';
tableNames{nEvals+3 +2}  = 'idx';
tableNames{nEvals+3 +3}  = 'Solver';
BIG_TABLE  = array2table(zeros( nConditions*nCases*nSolvers, nEvals+3+3 ), ...
  'VariableNames', tableNames);
count_keep = 0;

BIG_TABLE.SNR    = string(BIG_TABLE.SNR);
BIG_TABLE.Solver = string(BIG_TABLE.Solver);

%% MAIN LOOP
count_solv = 0;
count_solv_incr = 1/(nConditions*nCases*nSolvers);
WB_solv = waitbar(count_solv,'Sorting as a single table','Name','Collecting evaluations',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(WB_solv,'canceling',0);
for solverIDX = 1:nSolvers
  currSolver = SolveNames{solverIDX};
  for condition = 1:nConditions
    cd([ dataPath, '/', ConditionFolder{condition} ])
    % load the checklist and parameters
    load("evaluation.mat")
    % waitbar management
    waitbar(count_solv, WB_solv,...
      ['Solver : ',currSolver, '; SNR = ',num2str(info.SNRvals(condition)) ])
    count_case = 0;
    count_case_incr = 1/(nCases);
    %
    table_tmp = evaluation.(currSolver);
    table_tmp.SNR    = repmat( {num2str(info.SNRvals(condition))}, nCases, 1 );
    table_tmp.idx    = (1:nCases)';
    table_tmp.Solver = repmat( {currSolver}, nCases, 1 );
    %
    BIG_TABLE( count_keep*nCases + (1:nCases), : ) = table_tmp;
    count_keep = count_keep + 1;
    %
    if getappdata(WB_solv,'canceling')
      delete(WB_solv)
      break
    end
  end
  if getappdata(WB_solv,'canceling')
    delete(WB_solv)
    break
  end
end
delete(WB_solv)

% write final result
cd([originalPath, '/stats'])

writetable(TABLE,['EvalMetrics_' BaseName,'.xlsx']);