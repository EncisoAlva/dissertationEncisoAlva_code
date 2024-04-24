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
function RES = Protocol03_revised( meta, result, info )
RES = [];

idx    = 1:size(meta.Leadfield, 2);
meta.idxDips  = idx(( vecnorm( meta.Gridloc - meta.TrueCent, 2, 2 ) < 3*meta.kappa ));
meta.nDips    = length( meta.idxDips );

tmp = meta.idxDips*(3-1) + [1,2,3]';
meta.idxDips3 = sort(tmp(:));

flagSctr = false;
%debugCount = 0;
while ~flagSctr
  RES.SCtrIdx = randsample( meta.idxDips, 1);
  RES.SCtr    = meta.Gridloc(RES.SCtrIdx,:);
  if norm( RES.SCtr - meta.TrueCent ) < meta.kappa
    flagSctr = true;
  end
  %debugCount = debugCount + 1;
end
%disp(debugCount)
RES.tau  = rand(1,1)*(0.050-meta.kappa) + meta.kappa;
RES.idxS = idx(( vecnorm( meta.Gridloc - RES.SCtr, 2, 2 ) < 3*RES.tau ));

orient = randn(3,1);
orient = orient / norm(orient);

time   = linspace(0,1,30);
f_t    = ones(1, length(time));
weight = exp( -vecnorm( meta.Gridloc(meta.idxDips,:) - meta.TrueCent, 2, 2 ).^2/(2*(meta.kappa)^2) );

Yclean = meta.Leadfield(:,meta.idxDips3) * kron( orient, weight )*f_t;

PowY   = mean( Yclean.^2, 2 );

if isinf(meta.SNR)
  noise = zeros( size(Yclean) );
else
  noise  = normrnd(0, 1, size(Yclean) );
end

Y_OG   = Yclean + 10^(-meta.SNR/10) * diag(PowY) * noise;
Ywhite = lowpass(Y_OG, 2, 30); % post-whitening


RES.time   = time;
RES.Y_OG   = mean(Y_OG,2);
RES.Y      = mean(Ywhite,2);
RES.noise  = mean(noise,2);
RES.Yclean = mean(Yclean,2);

RES.vecY_OG   = (Y_OG);
RES.vecY      = (Ywhite);
RES.vecnoise  = (noise);
RES.vecYclean = (Yclean);

RES.TrueBarycent = weight' * meta.Gridloc(meta.idxDips,:) / sum(weight);

end