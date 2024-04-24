function collector(info)
% This function searches all the cases within my database, performs
% electrical source imaging with all the listed methods (which are tuned
% with the provided criteria), and then computes all the performance
% metrics over said reconstructions.
%
% This function is temporary. It was constructed for reproducibility and
% fast interaction with my database.
%
%-------------------------------------------------------------------------
% Author: Julio Cesar Enciso-Alva (2024)
%         juliocesar.encisoalva@mavs.uta.edu
%

BaseName   = info.BaseName;
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
cd([basePath, '/stats'])

writetable(BIG_TABLE,['EvalMetrics_' BaseName,'.xlsx']);

% back to the scripts folder
cd(originalPath)

end