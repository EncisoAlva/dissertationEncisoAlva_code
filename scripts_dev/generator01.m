% Creates a number of trials of synthetic data consisting on 
%  > Constrained sources at brain cortex
%  > 1 active dipole, random location
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
OGforward = 'asa_10_10';
BaseName  = 'test01';
nTrials   = 1000;

%% SETUP
originalPath = pwd;
addpath(pwd)
cd('..')
basePath  = pwd;
mkdir(BaseName)
cd(BaseName)

TotTrials = nTrials*6;
counter = 0;

%% PARAMETERS

% forward model
load([ basePath,'\og_data\', OGforward, '.mat']);

% metadata
meta = [];
meta.nDips = 1;

G = forward.Gain;
meta.nGridDips = size(G,2)/3;
meta.nChans    = size(G,1);
meta.Leadfield = zeros(meta.nChans, meta.nGridDips);
for i = 1:meta.nGridDips
  meta.Leadfield(:,i) = G( :, (3*(i-1)+(1:3)) ) * forward.GridOrient(i,:)';
end
meta.Gridloc   = forward.GridLoc;

save("metadata","meta");

%% MAIN LOOP
f = waitbar(counter,'','Name','Creating synthetic data',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(f,'canceling',0);
ResultPath = pwd;
for SNRi = [Inf, 40, 35, 30, 25, 20, 15, 10, 5, 0]
  cd(ResultPath)
  mkdir([ 'SNR_', num2str(SNRi) ])
  cd([ 'SNR_', num2str(SNRi) ])
  for idxRand = 1:nTrials
    meta.SNR = SNRi;
    % select source randomly
    rng(idxRand)
    meta.idxDips  = randi([1, meta.nGridDips]);
    meta.TrueLocs = forward.GridLoc(meta.idxDips,:);
    %
    RES = TestCase01( meta );
    %
    FileName = ['file',num2str(idxRand,'%04d'),'.mat'];
    %
    result = [];
    result.meta = meta;
    result.data = RES;
    result.meta.Leadfield = [];
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

% nested waitbars
if false
% Create bar for outer loop
h1 = waitbar(0,'I waitbar Please wait...');
for i=1:100
      % Create bar for inner loop
      h2=waitbar(0,'J waitbar Please wait...');
      % Change position of second bar so the is not overlap
      pos_w1=get(h1,'position');
      pos_w2=[pos_w1(1) pos_w1(2)+pos_w1(4) pos_w1(3) pos_w1(4)];
      set(h2,'position',pos_w2,'doublebuffer','on')
      % Scroll through first waitbar
      for j=1:100
          waitbar(j/100,h2)
      end
     % Scroll through second waitbar
      waitbar(i/100,h1)
      % Close first waitbar, recreate in each iteration of outer loop
      close(h2)
end
% Close final waitbar
close(h1)
end