function evaluator(info)
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

% which test data will be loaded
BaseName   = info.BaseName;

%% SETUP
originalPath = pwd;
addpath(pwd)
cd('..')
basePath = pwd;

%% INFO ABOUT TEST CASES
cd(['./data/', BaseName])
dataPath = pwd;

load("metadata.mat",'meta')
load("metadata2.mat",'info')

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
    load("checklist.mat",  'checklist')
    load("params.mat",     'params')
    load("evaluation.mat", 'evaluation')
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
      caseFile = ['file',num2str(0,'%04d'),'.mat'];
      load(caseFile,'result');
      %
      if info.print_all
        params.(currSolver) = feval( [currSolver,'_tune_pretty'] , meta, info, result );
      else
        params.(currSolver) = feval( [currSolver,'_tune'] , meta, info, result );
      end
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
      load(caseFile,'result');
      % actual solution
      solution = feval(currSolver, meta, info, result, params.(currSolver));
      % debug figures
      if info.debugFigs && (caseIDX < 2 )
        % REFERENCE SOURCES 
        Junit = zeros(meta.nGridDips,1);
        Junit(result.data.idxShort) = result.data.normJshort;
        %Junit = solution.normJ/max(solution.normJ(:));
        figure()
        trisurf(meta.Cortex.Faces, ...
          meta.Cortex.Vertices(:,1), meta.Cortex.Vertices(:,2), meta.Cortex.Vertices(:,3), ...
          'FaceColor', [1,1,1]*153/255, ...
          'EdgeColor', ...
          'none', 'FaceAlpha', 1 )
        view([ 90  90]) % top
        camlight('headlight', 'infinite')
        material dull
        %
        hold on
        trisurf(meta.Cortex.Faces, ...
          meta.Cortex.Vertices(:,1), meta.Cortex.Vertices(:,2), meta.Cortex.Vertices(:,3), ...
          'FaceColor', 'interp', ...
          'FaceVertexCData', Junit, ...
          'EdgeColor', 'none', ...
          'FaceAlpha', 'interp', ...
          'FaceVertexAlphaData', 1*(Junit>0.05) )
        material dull
        %colormap("turbo")
        colormap("parula")
        clim([0,1])
        %
        grid off
        set(gca,'DataAspectRatio',[1 1 1])
        set(gca,'XColor', 'none','YColor','none','ZColor','none')
        set(gca, 'color', 'none');
        set(gcf,'color','w');
        set(gca,'LooseInset',get(gca,'TightInset'))
        %
        b = colorbar;
        b.Label.String = 'Unitless; range=[0,1]';
        %
        fig = gcf;
        fig.Units = 'inches';
        fig.OuterPosition = [0 0 3 3];
        exportgraphics(gcf,[currSolver, '_ReferenceExample.pdf'],'Resolution',600)
        %
        % 50% 
        Junit = 1*(Junit > 0.5 );
        %
        figure()
        trisurf(meta.Cortex.Faces, ...
          meta.Cortex.Vertices(:,1), meta.Cortex.Vertices(:,2), meta.Cortex.Vertices(:,3), ...
          'FaceColor', [1,1,1]*153/255, ...
          'EdgeColor', ...
          'none', 'FaceAlpha', 1 )
        view([ 90  90]) % top
        camlight('headlight', 'infinite')
        material dull
        %
        hold on
        trisurf(meta.Cortex.Faces, ...
          meta.Cortex.Vertices(:,1), meta.Cortex.Vertices(:,2), meta.Cortex.Vertices(:,3), ...
          'FaceColor', 'interp', ...
          'FaceVertexCData', Junit, ...
          'EdgeColor', 'none', ...
          'FaceAlpha', 'interp', ...
          'FaceVertexAlphaData', 1*(Junit>0.05) )
        material dull
        %colormap("turbo")
        colormap("parula")
        clim([0,1])
        %
        grid off
        set(gca,'DataAspectRatio',[1 1 1])
        set(gca,'XColor', 'none','YColor','none','ZColor','none')
        set(gca, 'color', 'none');
        set(gcf,'color','w');
        set(gca,'LooseInset',get(gca,'TightInset'))
        %
        b = colorbar;
        b.Label.String = 'Unitless; range=[0,1]';
        %
        fig = gcf;
        fig.Units = 'inches';
        fig.OuterPosition = [0 0 3 3];
        exportgraphics(gcf,[currSolver, '_Reference_50perc_Example.pdf'],'Resolution',600)
        % ESTIMATED SOURCES 
        Junit = solution.normJ/max(solution.normJ(:));
        figure()
        trisurf(meta.Cortex.Faces, ...
          meta.Cortex.Vertices(:,1), meta.Cortex.Vertices(:,2), meta.Cortex.Vertices(:,3), ...
          'FaceColor', [1,1,1]*153/255, ...
          'EdgeColor', ...
          'none', 'FaceAlpha', 1 )
        view([ 90  90]) % top
        camlight('headlight', 'infinite')
        material dull
        %
        hold on
        trisurf(meta.Cortex.Faces, ...
          meta.Cortex.Vertices(:,1), meta.Cortex.Vertices(:,2), meta.Cortex.Vertices(:,3), ...
          'FaceColor', 'interp', ...
          'FaceVertexCData', Junit, ...
          'EdgeColor', 'none', ...
          'FaceAlpha', 'interp', ...
          'FaceVertexAlphaData', 1*(Junit>0.05) )
        material dull
        %colormap("turbo")
        colormap("parula")
        clim([0,1])
        %
        grid off
        set(gca,'DataAspectRatio',[1 1 1])
        set(gca,'XColor', 'none','YColor','none','ZColor','none')
        set(gca, 'color', 'none');
        set(gcf,'color','w');
        set(gca,'LooseInset',get(gca,'TightInset'))
        %
        b = colorbar;
        b.Label.String = 'Unitless; range=[0,1]';
        %
        fig = gcf;
        fig.Units = 'inches';
        fig.OuterPosition = [0 0 3 3];
        exportgraphics(gcf,[currSolver, '_EstimateExample.pdf'],'Resolution',600)
      end
      %
      % EVALUATION METRICS
      for evalIDX = 1:nEvals
        currEval = EvalNames{evalIDX};
        evaluation.(currSolver).(currEval)(caseIDX) = ...
          feval( currEval, meta, info, result, solution );
      end 
      meta.nGridDips
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

% back to the scripts folder
cd(originalPath)

end