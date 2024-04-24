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

% original forward model
OGforward  = 'asa_10_10';
OGatlas    = 'DesikanKilliany';
OGanatomy  = 'icbm152anatomy';

BaseName   = 'SimProt02';
nTrials    = 1000;
SourceType = 'surface';

%% SETUP
originalPath = pwd;
addpath(pwd)
cd('..')
basePath  = pwd;
mkdir(BaseName)
cd(BaseName)

TotTrials = nTrials*8;
counter = 0;

%% PARAMETERS

% forward model
load([ basePath,'\og_data\', OGforward, '.mat']);
load([ basePath,'\og_data\', OGatlas,   '.mat']);
load([ basePath,'\og_data\', OGanatomy, '.mat']);

% metadata
meta = [];
meta.nDips = 1;

G = forward.Gain;
meta.nRegions  = size( atlas, 2 );
meta.nGridDips = size(G,2)/3;
meta.nChans    = size(G,1);
meta.Leadfield = zeros(meta.nChans, meta.nGridDips);
for i = 1:meta.nGridDips
  meta.Leadfield(:,i) = G( :, (3*(i-1)+(1:3)) ) * forward.GridOrient(i,:)';
end
meta.Gridloc   = forward.GridLoc;
meta.atlas     = atlas;

% 'SIZE' OF DIPOLES
meta.GridVolume = zeros(meta.nGridDips,1);
switch SourceType
  case 'surface'
    for i = 1:size(cortex.Faces,1)
      tri  = cortex.Vertices( cortex.Faces(i,:), :);
      tri0 = tri([2,3],:) - tri(1,:);
      aar  = norm(cross( tri0(1,:), tri0(2,:) ))/2;
      meta.GridVolume(cortex.Faces(i,:)) = meta.GridVolume(cortex.Faces(i,:)) + aar/3;
    end
    % unit is m^2, please convert to mm^2 later
  case 'volume'
    DT = delaunayTriangulation(meta.Gridloc);
    for i = 1:size(DT.ConnectivityList,1)
      tet  = meta.Gridloc( DT.ConnectivityList(i,:), :);
      tet0 = tet([2,3,4],:) - tet(1,:);
      vvol = abs( cross( tet0(1,:),tet0(2,:) )*tet0(3,:)' )/6;
      meta.GridVolume(DT.ConnectivityList(i,:)) = meta.GridVolume(DT.ConnectivityList(i,:)) + vvol/4;
    end
end

save("metadata","meta");

%% MAIN LOOP
f = waitbar(counter,'','Name','Creating synthetic data',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(f,'canceling',0);
ResultPath = pwd;
for SNRi = [Inf, 35, 30, 25, 20, 15, 10, 5]
  cd(ResultPath)
  mkdir([ 'SNR_', num2str(SNRi) ])
  cd([ 'SNR_', num2str(SNRi) ])
  meta.SNR = SNRi;
  for idxRand = 1:nTrials
    % select source randomly
    rng(idxRand)
    meta.idxRegion = randi([1, meta.nRegions]);
    meta.Region    = atlas(meta.idxRegion).Label;
    meta.idxCent   = randsample( atlas(meta.idxRegion).Vertices, 1);
    meta.TrueCent  = forward.GridLoc(meta.idxCent,:);
    %meta.kappa     = (randsample( [0.01,.02,0.03],1 ));
    meta.kappa     = max(rand(1,1)*0.5, 0.001);
    %
    RES = SimulationProtocol02( meta );
    %
    FileName = ['file',num2str(idxRand,'%04d'),'.mat'];
    %
    result = [];
    result.meta = meta;
    result.data = RES;
    %
    % PRE-PREPROCESSING
    %SIG_12 = sqrtm(cov( result.data.Y ));
    %Ywhite = result.data.Y / SIG_12;
    %result.data.Y_OG   = result.data.Y; % original Y
    %result.data.Y      = Ywhite;        % Y after whitening
    %result.data.SIG_12 = SIG_12;        % whitening matrix
    %
    % computation of depth
    %result.data.TrueBarycent = result.meta.TrueLocs;
    [~,idxHead] = min(vecnorm( result.data.TrueBarycent - head.Vertices, 2, 2 ));
    [~,idxCort] = min(vecnorm( head.Vertices(idxHead,:) - cortex.Vertices, 2, 2 ));
    result.meta.Depth = norm( result.data.TrueBarycent - cortex.Vertices(idxCort,:) );
    %figure()
    %trisurf(cortex.Faces, cortex.Vertices(:,1), cortex.Vertices(:,2), cortex.Vertices(:,3), 'FaceAlpha', 0)
    %hold on
    %scatter3(result.data.TrueBarycent(1), result.data.TrueBarycent(2), result.data.TrueBarycent(3), 500, 'red')
    %
    result.meta.Leadfield  = [];
    result.meta.atlas      = [];
    result.meta.GridVolume = [];
    save(FileName,"result");
    %
    counter = counter + 1/TotTrials;
    waitbar(counter,f,['SNR=',num2str(SNRi)])
    if getappdata(f,'canceling')
      break
    end
  end
  if getappdata(f,'canceling')
    break
  end
end
delete(f)

%
cd(originalPath)

%writematrix(result.data.Y, 'TrashData.txt')
%Q = result.data.Y;
%save('TrashData.mat','Q')