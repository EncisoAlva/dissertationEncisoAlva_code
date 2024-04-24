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
function RES = TestCase01( meta )
RES = [];

time       = linspace(0,1,15);
SignalDips = sin(time*pi);

Yclean = meta.Leadfield(:,meta.idxDips) * ones(meta.nDips) * SignalDips;
PowY   = vecnorm(Yclean,2,2);

if isinf(meta.SNR)
  noise = zeros( size(Yclean) );
else
  noise  = normrnd(0, 1, size(Yclean) );
end
Y = Yclean + 10^(-meta.SNR/10) * diag(PowY) * noise;

RES.time   = time;
RES.Y      = Y;
RES.noise  = noise;
RES.Yclean = Yclean;

end