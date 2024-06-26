% It is meant to search for all test cases, solve each one with all
% available solver, and then run all evaluation metrics. The results from
% evaluation metrics are condensed on a big file.
%
% Author: Julio C Enciso-Alva (2024)
%         juliocesar.encisoalva@mavs.uta.edu
%

% original forward model
BaseName  = 'SimProt03';
%BaseName  = 'test03';
type = 'volume';

solversDone = {};
evalsDone   = {};

%% SETUP
originalPath = pwd;
addpath(pwd)
cd('..')
basePath  = pwd;

%% LOAD SOLVERS
cd('functions')
lst = dir('*.m');

nSolvers   = length(lst);
SolveNames = cell(nSolvers,1);
for i = 1:nSolvers
  tmp = split(lst(i).name,'.');
  SolveNames{i} = tmp{1};
end
% remove methods already done, running some of them is EXPENSIVE in time
SolveNames( ismember( SolveNames, solversDone ) ) = [];
nSolvers = length(SolveNames);

functionsPath = cd('..');
addpath(functionsPath)

%% LOAD EVALUATION METRICS
cd('test_functions')
lst = dir('*.m');

nEvals    = length(lst);
EvalNames = cell(nEvals,1);
for i = 1:nEvals
  tmp = split(lst(i).name,'.');
  EvalNames{i} = tmp{1};
end

evaluationsPath = cd('..');
addpath(evaluationsPath)

%% TEST CASES
cd(BaseName)

dataPath = pwd;
load('metadata.mat');
meta.type=type;

% first 2 folders are . and ..
lst = dir();
nConditions = -2;
for i =1:length(lst)
  if lst(i).isdir
    nConditions = nConditions+1;
  end
end
ConditionFolder = cell(nConditions,1);
counter = 3;
for i = 1:nConditions
  ConditionFolder{i} = lst(counter).name;
  counter = counter+1;
end

nCases = 500; % hard-coded for now

reset_result = false;

%% MAIN LOOP: SOLVING
counter = 0;
counter_incr = 1/(nConditions*nCases*nSolvers);
f = waitbar(counter,'(Starting)','Name','Evaluating methods | Outer Loop',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(f,'canceling',0);

for condition = 1:nConditions
  cd(ConditionFolder{condition})
  for solverIDX = 1:nSolvers
    currSolver = SolveNames{solverIDX};
    % sability patch
    %fileDebug  = ['debug_',currSolver,'.mat'];
    %if isfile(fileDebug)
    %  %
    %else
    %  completed_cases = 0;
    %end
    %
    f2 = waitbar(0,'Initial Parameter Tuning','Name',['Inner loop for ',currSolver],...
      'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    setappdata(f2,'canceling',0);
    % 2 extra evaluations: runtime of solver, runtime of parameter tuning
    waitbar(counter,f,['Condition: SNR=',ConditionFolder{condition},...
          ', Solver: ',currSolver])
    meta.optimizeParams = true;
    buffer = [];
    for caseIDX = 1:nCases
      waitbar(counter,f)
      waitbar(caseIDX/nCases,f2)
      counter = counter + counter_incr;
      caseFile = ['file',num2str(caseIDX,'%04d'),'.mat'];
      load(caseFile);
      meta.constrained = false;
      solution   = feval(currSolver, meta, result, buffer);
      if caseIDX==1
        meta.optimizeParams = false;
        buffer = solution.buffer;
      end
      solution.J = diag(meta.ColumnNorm.^-1) * solution.J;
      solutionFile = ['solution',num2str(caseIDX,'%04d'),'.mat'];
      if isfile(solutionFile)
        load(solutionFile);
      else
        collective_solution = [];
      end
      collective_solution.(currSolver) = solution;
      save(solutionFile,'collective_solution');
      if getappdata(f,'canceling')
        break
      end
      if getappdata(f2,'canceling')
        break
      end
    end
    delete(f2)
    %
    if getappdata(f,'canceling')
      break
    end
  end
  %
  cd(dataPath)
  if getappdata(f,'canceling')
    break
  end
end
delete(f)

reset_eval = true;

%% MAIN LOOP: EVALUATION
counter = 0;
counter_incr = 1/(nConditions*nCases*nSolvers);
f = waitbar(counter,'(Starting)','Name','Computing performance metrics',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(f,'canceling',0);
% pre-compute SVD decomposition of leadfield
%[U,S,V] = svd(meta.Leadfield);
%meta.U = U;
%meta.S = diag(S);
%meta.V = V;

for condition = 1:nConditions
  cd(ConditionFolder{condition})
  for solverIDX = 1:nSolvers
    currSolver = SolveNames{solverIDX};
    % load evaluation metrics already computed
    filename = ['evaluation_',ConditionFolder{condition},'.mat'];
    if isfile(filename)
      load(filename);
    else
      evaluation = [];
    end
    %
    % 2 extra evaluations: runtime of solver, runtime of parameter tuning
    if (~isfield(evaluation, currSolver)) || reset_eval
      evaluation.(currSolver) = zeros(nCases,3+nEvals);
    end
    if ~isfield(evaluation, 'SolveNames')
      evaluation.SolveNames = SolveNames;
      evaluation.EvalNames  = EvalNames;
      evaluation.ExtraNames = {'Algorithm Time','Parameter Tuning Time'};
    end
    for caseIDX = 1:nCases
      waitbar(counter,f,['Condition: SNR=',ConditionFolder{condition},...
        ', Solver: ',currSolver])
      counter = counter + counter_incr;
      caseFile = ['file',num2str(caseIDX,'%04d'),'.mat'];
      load(caseFile);
      solutionFile = ['solution',num2str(caseIDX,'%04d'),'.mat'];
      if isfile(solutionFile)
        load(solutionFile);
      else
        disp("Solutions file was not loaded.")
        continue
      end
      solution = collective_solution.(currSolver);
      %
      nExtras = 4;
      for evalIDX = 1:nEvals
        currEval = EvalNames{evalIDX};
        evaluation.(currSolver)(caseIDX, nExtras+evalIDX) =...
          feval( currEval, meta, result, solution );
      end
      evaluation.(currSolver)(caseIDX, 1) = solution.algTime;
      evaluation.(currSolver)(caseIDX, 2) = solution.parTime;
      evaluation.(currSolver)(caseIDX, 3) = result.meta.Depth*1000;
      evaluation.(currSolver)(caseIDX, 4) = result.meta.kappa*1000;
      if getappdata(f,'canceling')
        break
      end
    end
    %
    save(filename,'evaluation');
    if getappdata(f,'canceling')
      break
    end
  end
  %
  cd(dataPath)
  if getappdata(f,'canceling')
    break
  end
end
delete(f)

cd(originalPath)