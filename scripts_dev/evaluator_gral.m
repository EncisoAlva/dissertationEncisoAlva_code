% It is meant to search for all test cases, solve each one with all
% available solver, and then run all evaluation metrics. The results from
% evaluation metrics are condensed on a big file.
%
% Author: Julio C Enciso-Alva (2024)
%         juliocesar.encisoalva@mavs.uta.edu
%

cd(fileparts(which(mfilename)));
digits(60)

% which test data will be loaded
BaseName   = 'protocol_04_square500';
hard_reset = false;

%% SETUP
originalPath = pwd;
addpath(pwd)
cd('..')
basePath = pwd;

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
addpath('solvers', 'solvers_aux')

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
addpath('metrics', 'metrics_aux')

%% MAIN LOOP
count_solv = 0;
count_solv_incr = 1/(nConditions*nCases*nSolvers);
WB_solv = waitbar(count_solv,'Initializing','Name','Evaluating solvers',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(WB_solv,'canceling',0);
for solverIDX = 1:nSolvers
  currSolver = SolveNames{solverIDX};
  for condition = 1:nConditions
    cd([ dataPath, '/', ConditionFolder{condition} ])
    % load the checklist and parameters
    load("checklist.mat")
    load("params.mat")
    load("evaluation.mat")
    if ~isfield(params, currSolver)
      params.(currSolver) = [];
      save("params","params", '-v7.3');
    end
    if ~isfield(checklist, currSolver)
      checklist.(currSolver) = false(nCases,1);
      save("checklist","checklist")
    end
    if ~isfield(checklist.tuned, currSolver)
      checklist.tuned.(currSolver) = false;
      save("checklist","checklist")
    end
    if ~isfield(evaluation, currSolver)
      evaluation.(currSolver) = array2table(NaN(nCases, size(EvalNames,1)), ...
        'VariableNames',EvalNames);
      save("evaluation","evaluation")
    end
    % waitbar management
    waitbar(count_solv, WB_solv,...
      ['Solver : ',currSolver, '; SNR = ',num2str(info.SNRvals(condition)) ])
    count_case = 0;
    count_case_incr = 1/(nCases);
    WB_case = waitbar(0,'Initial Parameter Tuning','Name',...
      ['(Inner loop) ', currSolver, '; SNR = ',num2str(info.SNRvals(condition)) ],...
      'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    setappdata(WB_case,'canceling',0);
    %
    % perform parameter tuning ONCE and inform whenever it is complete
    if ~checklist.tuned.(currSolver)
      % select one or many cases at random
      caseFile = ['file',num2str(1,'%04d'),'.mat'];
      load(caseFile);
      %
      params.(currSolver) = feval( [currSolver,'_tune'] , meta, info, result );
      checklist.tuned.(currSolver) = true;
      save("params","params", '-v7.3');
      save("checklist","checklist");
    end
    waitbar(count_case, WB_case, 'Computing and Evaluating Solutions')
    for caseIDX = 1:nCases
      % waitbar update
      count_solv = count_solv + count_solv_incr;
      count_case = count_case + count_case_incr;
      waitbar(count_solv, WB_solv )
      waitbar(count_case, WB_case )
      % avoid reults already accounted for
      if checklist.(currSolver)(caseIDX)
        continue
      end
      % load individual data
      caseFile = ['file',num2str(caseIDX,'%04d'),'.mat'];
      load(caseFile);
      % actual solution
      solution = feval(currSolver, meta, info, result, params.(currSolver));
      % debug figures
      if info.debugFigs
        figure()
        trisurf(meta.Cortex.Faces, ...
          meta.Cortex.Vertices(:,1), meta.Cortex.Vertices(:,2), meta.Cortex.Vertices(:,3), 'FaceAlpha', 0)
        hold on
        scatter3(meta.Gridloc(:,1), meta.Gridloc(:,2), meta.Gridloc(:,3), ...
          40, solution.normJ*120,'filled')
        colormap("parula")
        scatter3( result.data.TrueCent(1), result.data.TrueCent(2), result.data.TrueCent(3), ...
          200, 'red','filled')
        title("Magnitude of true sources")
        b = colorbar;
        b.Label.String = 'Unitless; range=[0,120]';
      end
      %
      % EVALUATION METRICS
      for evalIDX = 1:nEvals
        currEval = EvalNames{evalIDX};
        evaluation.(currSolver).(currEval)(caseIDX) = ...
          feval( currEval, meta, info, result, solution );
      end
      evaluation.(currSolver).Run_Time(caseIDX) = ...
        solution.algTime;
      %evaluation.(currSolver).ParTune_Time(caseIDX) = ...
      %  solution.parTime;
      evaluation.(currSolver).Depth(caseIDX) = ...
        meta.GridDepth(result.idxCent);
      evaluation.(currSolver).Kappa(caseIDX) = ...
        result.kappa;
      %
      % track interrumpted executions
      checklist.(currSolver)(caseIDX) = true;
      %save("params","params", '-v7.3');
      save("checklist","checklist")
      save("evaluation","evaluation")
      if getappdata(WB_case,'canceling') || getappdata(WB_solv,'canceling')
        delete(WB_case);
        delete(WB_solv)
        break
      end
    end
    delete(WB_case)
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