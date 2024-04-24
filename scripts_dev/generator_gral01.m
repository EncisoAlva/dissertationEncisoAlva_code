% Creates a number of trials of synthetic data consisting on 
%  > Constrained sources at brain cortex
%  > 1 active region, only 1 dipole on that region
%  > Prescribed SNR, ranging from 0 to 40 in increments of 5, and also Inf
% Trials are saved as: 1 folder per SNR, 1 file per trial
%   > Folder (single SNR)
%   | > metadata, including leadfield
%   | > Trial 0001, only data
%   | > ...
%   | > Trial 1000
%
% Author: Julio C Enciso-Alva (2023)
%         juliocesar.encisoalva@mavs.uta.edu
%

%% PARAMETERS
info = [];

% forward model
info.OGforward  = 'asa_10_10_srf_BEM';
info.OGanatomy  = 'icbm152anatomy';

info.BaseName   = 'protocol_03_test';
info.SourceType = 'surface';

info.nTrials    = 50;
info.SNRvals    = [Inf, 30, 10, -10];

info.debugFigs  = false;

info.ProtocolFun   = '';
info.SourceProfile = 'square';

info.maxDepth  = 40; % unit: mm
info.maxKappa  = 30; % unit: mm
info.minKappa  =  1; % unit: mm

% for vol:  kap = 30.9 mm  ->  A = 30 cm^2
% for srf:  kap = 30.9 mm  ->  A = 30 cm^2

%% SETUP

scriptsPath = pwd;
addpath(pwd)
cd('..')
basePath  = pwd;
cd('./data')
mkdir(BaseName)
cd(BaseName)

TotTrials = nTrials*length(SNRvals);

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
meta.Type    = SourceType;
switch SourceType
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
switch SourceType
  case 'surface'
    ElementSize = zeros(size(meta.Cortex.Faces,1),1);
    for i = 1:size(meta.Cortex.Faces,1)
      tri  = meta.Cortex.Vertices( meta.Cortex.Faces(i,:), :);
      tri0 = tri([2,3],:) - tri(1,:);
      aar  = norm(cross( tri0(1,:), tri0(2,:) ))/2;
      meta.GridVolume(meta.Cortex.Faces(i,:)) = meta.GridVolume(meta.Cortex.Faces(i,:)) + aar/3;
      ElementSize(i) = aar;
    end
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
end

if debugFigs
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

% redundant informatoin
meta.info = info;

save("metadata","meta");

%% MAIN LOOP
counter = 0;
f = waitbar(counter,'','Name','Creating synthetic data',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(f,'canceling',0);
ResultPath = pwd;
for SNRi = SNRvals
  cd(ResultPath)
  mkdir([ 'SNR_', num2str(SNRi) ])
  cd([ 'SNR_', num2str(SNRi) ])
  result = [];
  result.SNR = SNRi;
  for idxRand = 1:nTrials
    % select source randomly
    rng(idxRand)
    found_source = false;
    while ~found_source
      tmp = randi([1 meta.nGridDips]);
      disp(tmp)
      if meta.GridDepth(tmp) < maxDepth
        found_source = true;
        result.idxCent = tmp;
      end
    end
    result.TrueCent  = meta.Gridloc(result.idxCent,:);
    result.kappa     = minKappa + rand(1,1)*(maxKappa-minKappa);
    %
    RES = SimulationProtocol03( meta, result );
    %
    FileName = ['file',num2str(idxRand,'%04d'),'.mat'];
    %
    result = [];
    result.data = RES;
    %
    % computation of depth
    [~,idxHead] = min(vecnorm( result.data.TrueBarycent - meta.Scalp.Vertices, 2, 2 ));
    [~,idxCort] = min(vecnorm( meta.Scalp.Vertices(idxHead,:) - meta.Cortex.Vertices, 2, 2 ));
    result.meta.Depth = norm( result.data.TrueBarycent - meta.Cortex.Vertices(idxCort,:) );
    %
    save(FileName,"result");
    %
    counter = counter + 1/TotTrials;
    waitbar(counter,f,['SNR=',num2str(SNRi)])
    if getappdata(f,'canceling')
      break
    end
  end
  % create empty files for solvers' parameters and checklists
  params    = [];
  checklist = [];
  save("params","params");
  save("checklist","checklist")
  %
  if getappdata(f,'canceling')
    break
  end
end
delete(f)

%
cd(scriptsPath)

%writematrix(result.data.Y, 'TrashData.txt')
%Q = result.data.Y;
%save('TrashData.mat','Q')