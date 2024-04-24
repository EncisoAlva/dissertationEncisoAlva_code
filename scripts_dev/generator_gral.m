% Creates a number of trials of synthetic data
%   > Folder (single SNR)
%   | > metadata, including leadfield
%   | > Trial 0001, only data
%   | > ...
%   | > Trial 1000
%
% Author: Julio C Enciso-Alva (2024)
%         juliocesar.encisoalva@mavs.uta.edu
%

cd(fileparts(which(mfilename)));
digits(60)

counter = 0;
f = waitbar(counter,'Computing general parameters, including SVD of leadfield matrix',...
  'Name','Computing synthetic data',...
  'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(f,'canceling',0);

%% PARAMETERS
info = [];

% forward model
info.OGforward  = 'asa_10_10_srf_BEM';
info.OGanatomy  = 'icbm152anatomy';

info.BaseName   = 'protocol_04_square500';
info.SourceType = 'surface';

info.nTrials    = 500;
info.SNRvals    = [Inf, 30, 20, 10, 0, -10];

info.ProtocolFun   = 'Protocol04';
info.SourceProfile = 'square';

info.maxDepth  = 25; % unit: mm
info.maxKappa  = 10*sqrt(10/pi); % unit: mm
info.minKappa  = 10*sqrt(10/pi); % unit: mm

% for vol:  kap = 30.9 mm  ->  A = 30 cm^2
% for srf:  kap = 30.9 mm  ->  A = 30 cm^2

info.debugFigs  = false;

info.debugCent  = true;
info.debugCoord = [47.353, 18.555, 113.019];

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
load([ basePath,'\anat_ref\', info.OGforward, '.mat']);
load([ basePath,'\anat_ref\', info.OGanatomy, '.mat']);

% metadata
meta = [];

% adding cortex and scalp triangulations
meta.Cortex.Vertices = cortex.Vertices*1000;
meta.Cortex.Faces    = cortex.Faces;
meta.Cortex.Normals  = cortex.VertNormals;
meta.Scalp.Vertices  = head.Vertices*1000;
meta.Scalp.Faces     = head.Faces;
meta.Scalp.Normals   = head.VertNormals;

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
end

% some figres related to dipoles' depth and area/volume
if info.debugFigs
  figure()
  trisurf(meta.Cortex.Faces, ...
    meta.Cortex.Vertices(:,1), meta.Cortex.Vertices(:,2), meta.Cortex.Vertices(:,3), 'FaceAlpha', 0)
  hold on
  scatter3(meta.Gridloc(:,1), meta.Gridloc(:,2), meta.Gridloc(:,3), ...
    40, meta.GridDepth,'filled')
  colormap("parula")
  title('Depth of dipoles: distance from scalp')
  a=colorbar;
  a.Label.String = 'Distance to scalp [mm]';
  %
  figure()
  clf
  trisurf(meta.Cortex.Faces, ...
    meta.Cortex.Vertices(:,1), meta.Cortex.Vertices(:,2), meta.Cortex.Vertices(:,3), 'FaceAlpha', 0)
  hold on
  scatter3(meta.Gridloc(:,1), meta.Gridloc(:,2), meta.Gridloc(:,3), ...
    40, meta.GridVolume,'filled')
  colormap("parula")
  b=colorbar;
  switch meta.Type
    case 'surface'
      title("'Area' of dipoles: area of neighboring triangle elements")
      b.Label.String = 'Distance to scalp [mm^2]';
    case 'volume'
      title("'Volume' of dipoles: volume of neighboring tetahedra elements")
      b.Label.String = 'Distance to scalp [mm^3]';
  end
end

% for geodesic distance
switch info.SourceType
  case 'surface'
    meta.DipEdges = edges_mine( meta.Cortex );
  case 'volume'
    meta.DipEdges = edges_mine( meta.DT );
end
meta.asGraph = graph( table(meta.DipEdges, 'VariableNames',{'EndNodes'}) );
minDist = zeros(meta.nGridDips,1);
for ii = 1:meta.nGridDips
  minDist(ii) = min(vecnorm( ...
    meta.Gridloc([(1:(ii-1)),((ii+1):meta.nGridDips)],:) - meta.Gridloc(ii,:), 2, 2 ));
end
meta.minDist = minDist;

% column-normalization
meta.ColumnNorm = vecnorm(meta.LeadfieldOG,2,1);
meta.Leadfield  = meta.LeadfieldOG;
for q = 1:size(meta.Leadfield,2)
  meta.Leadfield(:,q) = meta.Leadfield(:,q)/meta.ColumnNorm(q);
end

% SVD decomposition of leadfield matrix
[U,S,V] = svd(meta.Leadfield,'vector');
meta.U  = U;
meta.S  = S;
meta.V  = V;

% redundancy
meta.info = info;
save("metadata","meta");
save("metadata2","info");

%% MAIN LOOP
ResultPath = pwd;
for SNRi = info.SNRvals
  cd(ResultPath)
  mkdir([ 'SNR_', num2str(SNRi) ])
  cd([ 'SNR_', num2str(SNRi) ])
  result = [];
  result.SNR = SNRi;
  for idxRand = 1:info.nTrials
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
    % bypass selection of a random center to a known point
    if info.debugCent
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
    if info.debugFigs
      figure()
      clf
      trisurf(meta.Cortex.Faces, ...
        meta.Cortex.Vertices(:,1), meta.Cortex.Vertices(:,2), meta.Cortex.Vertices(:,3), 'FaceAlpha', 0)
      hold on
      scatter3(result.IntendedCent(1), result.IntendedCent(2), result.IntendedCent(3), ...
        300, 'red','filled')
    end
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