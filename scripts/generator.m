function generator(info)
% This function creates a set of test cases based on the sepcifications
% included in the structure info.
%
% This function is temporary. It was constructed for reproducibility and
% fast interaction with my database.
%
%-------------------------------------------------------------------------
% Author: Julio Cesar Enciso-Alva (2024)
%         juliocesar.encisoalva@mavs.uta.edu
%

counter = 0;
f = waitbar(counter,'Computing general parameters, including SVD of leadfield matrix',...
  'Name','Computing synthetic data',...
  'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(f,'canceling',0);

%% SETUP

scriptsPath = pwd;
addpath(pwd)
cd('..')
basePath  = pwd;
cd('./data')
mkdir(info.BaseName)
cd(info.BaseName)

TotTrials = info.nTrials*length(info.SNRvals);

%% LOAD DATA

% forward model
load([ basePath,'\anat_ref\', info.OGforward, '.mat'],'forward');
load([ basePath,'\anat_ref\', info.OGanatomy, '.mat'],'cortex','head');

% metadata
meta = [];

% adding cortex and scalp triangulations
meta.Cortex.Vertices = cortex.Vertices*1000;
meta.Cortex.Faces    = cortex.Faces;
meta.Cortex.Normals  = cortex.VertNormals;
meta.Scalp.Vertices  = head.Vertices*1000;
meta.Scalp.Faces     = head.Faces;
meta.Scalp.Normals   = head.VertNormals;

% remove removed channels (absent by default, didn't noticed before)
forward.Gain(any(isnan(forward.Gain), 2), :) = [];

% leadfield matrix and general info
meta.nChans  = size(forward.Gain,1);
meta.Gridloc = forward.GridLoc*1000; % standard unit is m, changed to mm
meta.nGridDips = size(forward.Gain,2)/3;
meta.Type    = info.SourceType;
switch info.SourceType
  case 'surface'
    meta.LeadfieldOG = zeros(meta.nChans, meta.nGridDips);
    for i = 1:meta.nGridDips
      meta.LeadfieldOG(:,i) = forward.Gain( :, (3*(i-1)+(1:3)) ) * forward.GridOrient(i,:)';
    end
  case 'volume'
    meta.LeadfieldOG = forward.Gain;
end

% depth of dipoles: distance to scalp
meta.GridDepth = zeros(meta.nGridDips,1);
for ii = 1:meta.nGridDips
  meta.GridDepth(ii) = min(vecnorm( meta.Gridloc(ii,:) - meta.Scalp.Vertices, 2, 2 ));
end

% surface/volume of dipoles: sum of surf/vol of elements
meta.GridVolume = zeros(meta.nGridDips,1);
switch info.SourceType
  case 'surface'
    ElementSize = zeros(size(meta.Cortex.Faces,1),1);
    for i = 1:size(meta.Cortex.Faces,1)
      tri  = meta.Cortex.Vertices( meta.Cortex.Faces(i,:), :);
      tri0 = tri([2,3],:) - tri(1,:);
      aar  = norm(cross( tri0(1,:), tri0(2,:) ))/2;
      meta.GridVolume(meta.Cortex.Faces(i,:)) = meta.GridVolume(meta.Cortex.Faces(i,:)) + aar/3;
      ElementSize(i) = aar;
    end
    meta.Elements    = meta.Cortex.Faces;
    meta.ElementSize = ElementSize;
  case 'volume'
    DT = delaunayTriangulation(meta.Gridloc);
    nElements   = size(DT.ConnectivityList,1);
    % find 'faulty' elemts that are outside the volume
    ElementCent = zeros(nElements,3);
    for ii = 1:nElements
      ElementCent(ii,:) = mean( meta.Gridloc(DT(ii,:)', :), 1);
    end
    insideElement = inpolyhedron(meta.Cortex.Faces, meta.Cortex.Vertices, ElementCent);
    cleanConnectivityList = DT.ConnectivityList(insideElement,:);
    %
    ElementSize = zeros(size(cleanConnectivityList,1),1);
    for i = 1:size(cleanConnectivityList,1)
      tet  = meta.Gridloc( cleanConnectivityList(i,:), :);
      tet0 = tet([2,3,4],:) - tet(1,:);
      vvol = abs( cross( tet0(1,:),tet0(2,:) )*tet0(3,:)' )/6;
      meta.GridVolume(cleanConnectivityList(i,:)) = meta.GridVolume(cleanConnectivityList(i,:)) + vvol/4;
      ElementSize(i) = vvol;
    end
    meta.Elements    = cleanConnectivityList;
    meta.ElementSize = ElementSize;
    meta.DT = DT;
end

% some figres related to dipoles' depth and area/volume
if info.debugFigs
  figure()
  trisurf(meta.Cortex.Faces, ...
    meta.Cortex.Vertices(:,1), meta.Cortex.Vertices(:,2), meta.Cortex.Vertices(:,3), ...
    'FaceColor', 'interp', ...
    'FaceVertexCData', meta.GridDepth, ...
    'EdgeColor', 'none' )
  view([ 90  90]) % top
  camlight('headlight', 'infinite')
  material dull
  %
  set(gca,'DataAspectRatio',[1 1 1])
  %set(gca,'XColor', 'none','YColor','none','ZColor','none')
  set(gca,'XColor', 'black','YColor','black','ZColor','black')
  set(gca, 'color', 'none');
  %grid off
  %title('Depth of dipoles: distance from scalp')
  a=colorbar;
  a.Label.String = 'Distance to scalp [mm]';
  %
  xlabel('[mm]')
  ylabel('[mm]')
  zlabel('[mm]')
  %view([  0   0]) % right
  %view([ 90   0]) % front
  %view([180   0]) % left
  %view([270   0]) % back
  %view([ 90  90]) % top
  %
  figure()
  trisurf(meta.Cortex.Faces, ...
    meta.Cortex.Vertices(:,1), meta.Cortex.Vertices(:,2), meta.Cortex.Vertices(:,3), ...
    'FaceColor', 'interp', ...
    'FaceVertexCData', meta.GridVolume, ...
    'EdgeColor', 'none' )
  view([ 90  90]) % top
  camlight('headlight', 'infinite')
  material dull
  %
  set(gca,'DataAspectRatio',[1 1 1])
  %set(gca,'XColor', 'none','YColor','none','ZColor','none')
  %set(gca,'XColor', 'black','YColor','black','ZColor','black')
  set(gca, 'color', 'none');
  %
  b=colorbar;
  switch meta.Type
    case 'surface'
      %title("Area of neighboring elements")
      b.Label.String = "Area of neighboring elements [mm^2]";
    case 'volume'
      %title("'Volume' of dipoles: volume of neighboring tetahedra elements")
      b.Label.String = "Voluma of neighboring elements [mm^3]";
  end
  %
  xlabel('[mm]')
  ylabel('[mm]')
  zlabel('[mm]')
  %view([  0   0]) % right
  %view([ 90   0]) % front
  %view([180   0]) % left
  %view([270   0]) % back
  %view([ 90  90]) % top
end

% for geodesic distance
switch info.SourceType
  case 'surface'
    meta.DipEdges = edges_mine( meta.Cortex );
    weights = zeros( size(meta.DipEdges,1), 1);
    for i = 1:size(meta.DipEdges,1)
      weights(i) = norm( meta.Cortex.Vertices(meta.DipEdges(i,1),:) - meta.Cortex.Vertices(meta.DipEdges(i,2),:), 2 );
    end
    %meta.asGraph = graph( table(meta.DipEdges, 'VariableNames',{'EndNodes'}) );
    tmp_table = table(meta.DipEdges, 'VariableNames',{'EndNodes'});
    tmp_table.Weight = weights;
    meta.asGraph = graph( tmp_table );
  case 'volume'
    %meta.DipEdges = edges_mine( meta.DT );
    %weights = zeros( size(meta.DipEdges,1), 1);
    %for i = 1:size(meta.DipEdges,1)
    %  weights(i) = norm( meta.DT.Vertices(meta.DipEdges(i,1),:) - meta.DT.Vertices(meta.DipEdges(i,2),:), 2 );
    %end
end
%minDist = zeros(meta.nGridDips,1);
%for ii = 1:meta.nGridDips
%  minDist(ii) = min(vecnorm( ...
%    meta.Gridloc([(1:(ii-1)),((ii+1):meta.nGridDips)],:) - meta.Gridloc(ii,:), 2, 2 ));
%end
%meta.minDist = minDist;

% re-referencing
meta.LeadfieldAvg = meta.LeadfieldOG - mean( meta.LeadfieldOG, 1 );
if info.debugFigs
  figure()
  tiledlayout(1,2)
  nexttile
  histogram(sum(meta.LeadfieldOG, 1))
  title('1 * G(:,i) ')
  nexttile
  histogram(sum(meta.LeadfieldAvg, 1))
  title('1 * G(:,i) ')
end

% column-normalization
meta.ColumnNorm       = vecnorm(meta.LeadfieldAvg,2,1);
meta.Leadfield        = meta.LeadfieldAvg;
meta.LeadfieldColNorm = meta.LeadfieldAvg;
for q = 1:size(meta.Leadfield,2)
  meta.LeadfieldColNorm(:,q) = meta.LeadfieldColNorm(:,q)/meta.ColumnNorm(q);
end

% SVD decomposition of leadfield matrix
[U,S,V] = svd(meta.Leadfield,'vector');
meta.U  = U;
meta.S  = S;
meta.V  = V;

% redundancy
meta.info = info;
save("metadata","meta", "-v7.3");
save("metadata2","info");

%% MAIN LOOP
ResultPath = pwd;
for SNRi = info.SNRvals
  cd(ResultPath)
  mkdir([ 'SNR_', num2str(SNRi) ])
  cd([ 'SNR_', num2str(SNRi) ])
  result = [];
  result.SNR = SNRi;
  for idxRand = 0:info.nTrials
    if idxRand == 0
      % curated points for parameter tuning
      [~,result.idxCent] = min(vecnorm( meta.Gridloc - info.debugCoord, 2, 2));
    else
      % select source randomly
      rng(idxRand)
      found_source = false;
      while ~found_source
        tmp = randi([1 meta.nGridDips]);
        if info.debugFigs
          disp(tmp)
        end
        if meta.GridDepth(tmp) < info.maxDepth
          found_source = true;
          result.idxCent = tmp;
        end
      end
    end
    % bypass selection of a random center to a known point
    if info.debugCent
      % curated points for parameter tuning
      [~,result.idxCent] = min(vecnorm( meta.Gridloc - info.debugCoord, 2, 2));
    end
    result.IntendedCent = meta.Gridloc(result.idxCent,:);
    tmp = randn(3,1);
    result.Orient   = tmp / norm(tmp);
    % extension is either (1) a given value, or (2) randomized within a range
    if info.minKappa == info.maxKappa
      result.kappa = info.minKappa;
    else
      result.kappa = info.minKappa + rand(1,1)*(info.maxKappa-info.minKappa);
    end
    %
    result.data = feval(info.ProtocolFun, meta, result, info);
    %
    FileName = ['file',num2str(idxRand,'%04d'),'.mat'];
    save(FileName,"result");
    %
    counter = counter + 1/TotTrials;
    waitbar(counter,f,['SNR=',num2str(SNRi)])
    if getappdata(f,'canceling')
      break
    end
  end
  % create empty files for solvers' parameters and checklists
  params     = [];
  evaluation = [];
  checklist  = []; 
  checklist.tuned = [];
  save("params","params");
  save("evaluation","evaluation")
  save("checklist","checklist")
  %
  if getappdata(f,'canceling')
    break
  end
end
delete(f)

% back to the scripts folder
cd(scriptsPath)

end