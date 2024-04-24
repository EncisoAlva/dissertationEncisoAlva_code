% This script creates a single trials of synthetic data according to 
% protocol in the Multiple Source Preallocation Paper:
%  > Constrained dipoles at brain cortex
%  > One single active dipole (given)
%  > Sample freq = 15 Hz [actually irrelevant]
%  > Sample window = 1 sec, from 0 to 1
%  > Total points: 15
%  > Signal: Sine wave with peak at t = 0.5 s
%  > Added noise on sensors with prescribed SNR (given)
%
% Author: Julio C Enciso-Alva (2023)
%         juliocesar.encisoalva@mavs.uta.edu
%
function RES = Protocol04( meta, result, info )

RES = [];

% optional: only consider sources with magnitude > 5%
% the maximal draw distance depend on the profile
switch info.SourceProfile
  case 'square'
    maxDist = result.kappa;
  case 'exp'
    maxDist = 2.31*result.kappa;
  case 'gauss'
    maxDist = 1.52*result.kappa;
  case 'circ'
    maxDist = result.kappa;
end
% prepare a short list of dipoles within the draw distance
idx = 1:size(meta.Leadfield, 2);
RES.idxShort  = idx( vecnorm( meta.Gridloc - result.IntendedCent, 2, 2 ) < maxDist*1.5 );
RES.nShort    = length( RES.idxShort );
RES.idxCentShort = find(result.idxCent == RES.idxShort,1);
if isempty(RES.idxCentShort)
  [~,RES.idxCentShort] = min(vecnorm( meta.Gridloc(RES.idxShort,:) - result.IntendedCent, 2, 2));
end

% miscellanea
GraphDist = distances(meta.asGraph, RES.idxShort, RES.idxShort) * max(meta.minDist(RES.idxShort));
switch info.SourceType
  case 'volume'
    tmp = (RES.idxShort-1)*3 + [1,2,3]';
    RES.idxShortG = sort(tmp(:));
  case 'surface'
    RES.idxShortG = RES.idxShort;
end

% J, shortened to non-zero dipoles
RES.time   = linspace(0,0,1);
RES.Jshort = zeros(RES.nShort,1);
switch info.SourceProfile
  case 'square'
    for ii = 1:RES.nShort
      switch info.SourceType
        case 'surface'
          % geodesic distance in cortex
          if GraphDist(ii,RES.idxCentShort) < result.kappa
            RES.Jshort(ii) = 1;
          end
        case 'volume'
          % euclidian distance in 3D space
          idx = RES.idxShort(ii); 
          if vecnorm( meta.Gridloc(idx,:) - result.IntendedCent, 2, 2 ) < result.kappa
            RES.Jshort(ii) = 1;
          end
      end
    end
  case 'exp'
    switch info.SourceType
      case 'surface'
        RES.Jshort = exp(- GraphDist(:,RES.idxCentShort) /result.kappa);
      case 'volume'
        RES.Jshort = exp(-vecnorm( meta.Gridloc(RES.idxShort,:) - result.IntendedCent, 2, 2 )/result.kappa);
    end
  case 'gauss'
    switch info.SourceType
      case 'surface'
        RES.Jshort = exp(-( GraphDist(:,RES.idxCentShort) ).^2/(2*(result.kappa^2)));
      case 'volume'
        RES.Jshort = exp(-vecnorm( meta.Gridloc(RES.idxShort,:) - result.IntendedCent, 2, 2 ).^2/(2*(result.kappa^2)));
    end
  case 'circ'
    switch info.SourceType
      case 'surface'
        RES.Jshort = ( 1 - ( GraphDist(:,RES.idxCentShort) /result.kappa ).^2 ).^(1/2);
      case 'volume'
        RES.Jshort = ( 1 - ...
          ( vecnorm( meta.Gridloc(RES.idxShort,:) - result.IntendedCent, 2, 2 ) /result.kappa ).^2 ).^(1/2);
    end
end

% the norm of J, also shortened to save space
switch info.SourceType
  case 'surface'
    RES.normJshort = abs( RES.Jshort );
  case 'volume'
    RES.normJshort = dip_norm( RES.Jshort );
end

% debug figures
if info.debugFigs
  figure()
  trisurf(meta.Cortex.Faces, ...
    meta.Cortex.Vertices(:,1), meta.Cortex.Vertices(:,2), meta.Cortex.Vertices(:,3), 'FaceAlpha', 0)
  hold on
  scatter3(meta.Gridloc(RES.idxShort,1), meta.Gridloc(RES.idxShort,2), meta.Gridloc(RES.idxShort,3), ...
    40, RES.normJshort*120,'filled')
  colormap("parula")
  scatter3( result.IntendedCent(1), result.IntendedCent(2), result.IntendedCent(3), ...
    200, 'red','filled')
  title("Magnitude of true sources")
  b = colorbar;
  b.Label.String = 'Unitless; range=[0,120]';
end

% Y, noiseless
switch info.SourceType
  case 'volume'
    RES.Yclean = meta.LeadfieldOG(:,RES.idxShortG) * kron( result.Orient, RES.Jshort );
  case 'surface'
    RES.Yclean = meta.LeadfieldOG(:,RES.idxShortG) * RES.Jshort;
end
%RES.varY = vecnorm( meta.LeadfieldOG, 2, 2 );
RES.varY = RES.Yclean.^2;

% adding noise to a prescribed SNR
if isinf(result.SNR)
  noise = zeros( size(RES.Yclean) );
else
  noise = normrnd(0, 1, size(RES.Yclean) );
end

% Y
RES.YOG = RES.Yclean + 10^(-result.SNR/10) * diag( RES.varY ) * noise;
RES.Y   = RES.YOG - mean(RES.YOG,1);

% true center of mass
RES.TrueCent = RES.Jshort' * meta.Gridloc(RES.idxShort,:) / sum(RES.Jshort);

end