function RES = TestCase03( meta )

RES = [];

time       = linspace(0,1,15);
SignalDips = sin(time*pi);

orient = randn(3,1);
orient = orient / norm(orient);

Yclean = meta.Leadfield(:, (3*(meta.idxDips-1) + (1:3)) ) * orient * ...
  ones(meta.nDips) * SignalDips;
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
RES.orient = orient;

end