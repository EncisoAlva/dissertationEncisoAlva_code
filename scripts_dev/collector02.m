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

% original forward model
BaseName  = 'SimProt02';
%BaseName  = 'test03';

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

nCases = 1000; % hard-coded for now

reset_eval = false;
TABLE = table;

%% MAIN LOOP: EVALUATION
counter = 0;
counter_incr = 1/(nConditions*nSolvers);
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
  tmpLabel = split( ConditionFolder{condition} ,'_') ;
  SNRlabel = tmpLabel{2};
  for solverIDX = 1:nSolvers
    currSolver = SolveNames{solverIDX};
    % load evaluation metrics already computed
    filename = ['evaluation_',ConditionFolder{condition},'.mat'];
    if ~isfile(filename)
      continue
    end
    load(filename)
    if ~isfield(evaluation, currSolver)
      continue
    end
    %
    % 2 extra evaluations: runtime of solver, runtime of parameter tuning
    % 4 extra data: solvername, SNR, depth of true barycenter, kappa
    data = zeros(nCases,nEvals +2 +4);
    data(:, 2+(1:(nEvals+4)) ) = evaluation.(currSolver);
    TABLEtmp = array2table(data, ...
      'VariableNames', cat(2,{'Solver', 'SNR','Depth','Kappa'},evaluation.ExtraNames,EvalNames') );
    TABLEtmp.Solver = repmat( convertCharsToStrings(currSolver),nCases,1);
    TABLEtmp.SNR    = repmat( convertCharsToStrings(SNRlabel),  nCases,1);
    if max( condition,solverIDX ) == 1
      TABLE = TABLEtmp;
    else
      TABLE = [TABLE; TABLEtmp];
    end
    %
    waitbar(counter,f,['Condition: SNR=',ConditionFolder{condition},...
      ', Solver: ',currSolver])
    counter = counter + counter_incr;
    if getappdata(f,'canceling')
      break
    end
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

writetable(TABLE,['EvalMetrics_' BaseName,'.xlsx']);

% LocErr = LocErr(:,[8,6,5,4,3,2,1,7]);
% HMSize = HMSize(:,[8,6,5,4,3,2,1,7]);
% 
% ConditionFolder2 = cell(8,1);
% ConditionFolder2{1} = ConditionFolder{8};
% ConditionFolder2{2} = ConditionFolder{6};
% ConditionFolder2{3} = ConditionFolder{5};
% ConditionFolder2{4} = ConditionFolder{4};
% ConditionFolder2{5} = ConditionFolder{3};
% ConditionFolder2{6} = ConditionFolder{2};
% ConditionFolder2{7} = ConditionFolder{1};
% ConditionFolder2{8} = ConditionFolder{7};
% 
% for i = 1:8
%   tmp = split(ConditionFolder2{i},'_');
%   ConditionFolder2{i} = tmp{2};
% end
% 
% figure()
% boxplot(100*LocErr)
% xlabel('SNR [dB]')
% ylabel('Localization Error [cm]')
% xticklabels(ConditionFolder2)
% grid on
% title('Tikhonov')
% 
% fig = gcf;
% fig.Color = [1,1,1];
% 
% %export_fig 'boxplot_LocErr' -pdf -r300 
% 
% figure()
% boxplot(HMSize)
% xlabel('SNR [dB]')
% ylabel('Half-Max Size [dipoles]')
% grid on
% xticklabels(ConditionFolder2)
% 
% fig = gcf;
% fig.Color = [1,1,1];
% 
% %export_fig 'boxplot_HMS' -pdf -r300 

cd(originalPath)