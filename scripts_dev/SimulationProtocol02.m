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
function RES = SimulationProtocol02( meta )
RES = [];

idx    = 1:size(meta.Leadfield, 2);
meta.idxDips = idx(( vecnorm( meta.Gridloc - meta.TrueCent, 2, 2 ) < 4*meta.kappa ));
meta.nDips   = length( meta.idxDips );

weight = exp( vecnorm( -meta.Gridloc(meta.idxDips,:) - meta.TrueCent, 2, 2 ).^2/(2*meta.kappa^2) );

time   = linspace(0,1,15);
f_t    = sin(time*pi);

Yclean = meta.Leadfield(:,meta.idxDips) * weight * f_t;
%PowY   = vecnorm(Yclean,2,2)^2;
PowY   = mean( Yclean.^2, 2 );

if isinf(meta.SNR)
  noise = zeros( size(Yclean) );
else
  noise  = normrnd(0, 1, size(Yclean) );
end

Y      = Yclean + 10^(-meta.SNR/10) * diag(PowY) * noise;
Ywhite = Yclean + noise; % post-whitening


RES.time   = time;
RES.Y_OG   = Y;
if isinf(meta.SNR)
  RES.SIG_12 = eye( length(PowY) );
else
  RES.SIG_12 = 10^(-meta.SNR/20) * diag(sqrt( PowY ));
end
RES.Y      = inv(RES.SIG_12) * Ywhite;
RES.noise  = noise;
RES.Yclean = Yclean;

RES.TrueBarycent = weight' * meta.Gridloc(meta.idxDips,:) / sum(weight);

end