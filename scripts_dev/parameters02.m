% Pre-computes some information from the matrices that are common to all
% trials, in order to save time (but spending space on disk)
%  > Column normalized leadfield matrices
%
% Author: Julio C Enciso-Alva (2023)
%         juliocesar.encisoalva@mavs.uta.edu
%

%TargetDirs = {'test02','test_vol01'};SimProt02
%TargetDirs = {'test03'};
TargetDirs = {'SimProt03'};
SourceType = {'volume'};

%% SETUP
originalPath = pwd;
addpath(pwd)
cd('..')
BasePath  = pwd;

counter = 0;
incr    = 1/(size(TargetDirs,2));

reset = false;

%% MAIN LOOP

f = waitbar(counter,'','Name','Pre-computing parameters',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(f,'canceling',0);
for i = 1:size(TargetDirs,2)
BaseName  = TargetDirs{i};
waitbar(counter,f,['Current directory : ',BaseName])
counter = counter + incr;
if getappdata(f,'canceling')
  break
end

cd(BaseName)
load("metadata.mat");

% COLUMN NORMALIZATION
if ~isfield(meta,'LeadfieldOG') || reset
  meta.ColumnNorm  = vecnorm(meta.Leadfield,2,1);
  meta.LeadfieldOG = meta.Leadfield;
  for q = 1:size(meta.Leadfield,2)
    meta.Leadfield(:,q) = meta.Leadfield(:,q)/meta.ColumnNorm(q);
  end
end

% PARAMETERS

% SVD decomposition of leadfield matrix
if ~isfield(meta,'U') || reset
  [U,S,V] = svd(meta.Leadfield,'vector');
  meta.U  = U;
  meta.S  = S;
  meta.V  = V;
end

% save and go back to the next
save('metadata','meta', '-v7.3');
cd(BasePath)
end
delete(f)

cd(originalPath)