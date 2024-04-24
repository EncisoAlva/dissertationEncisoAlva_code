% Pre-computes some information from the matrices that are common to all
% trials, in order to save time (but spending space on disk)
%  > Column normalized leadfield matrices
%
% Author: Julio C Enciso-Alva (2023)
%         juliocesar.encisoalva@mavs.uta.edu
%

TargetDirs = {'test01', 'test02', 'test_vol01'};

%% SETUP
originalPath = pwd;
addpath(pwd)
cd('..')
BasePath  = pwd;

counter = 0;
incr    = 1/size(TargetDirs,2);

%% MAIN LOOP

f = waitbar(counter,'','Name','Pre-computing parameters',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(f,'canceling',0);
for i = 1:size(TargetDirs,2)
BaseName  = TargetDirs{i};
waitbar(counter,f,['Current directory : ',BaseName])
counter = counter + incr;

cd(BaseName)
load("metadata.mat");

% PARAMETERS

% SVD decomposition of leadfield matrix
[U,S,V] = svd(meta.Leadfield);
meta.U  = U;
meta.S  = diag(S);
meta.V  = V;

% row-normalized
meta.ColNorm.Norm      = vecnorm(meta.Leadfield,2,1);
meta.ColNorm.Leadfield = zeros( size(meta.Leadfield) );
for col = 1:size(meta.ColNorm.Leadfield, 2)
  meta.ColNorm.Leadfield(:,col) = meta.Leadfield(:,col) / meta.ColNorm.Norm(col);
end
[U,S,V] = svd(meta.ColNorm.Leadfield);
meta.ColNorm.U  = U;
meta.ColNorm.S  = diag(S);
meta.ColNorm.V  = V;

% save and go back to the next
save("metadata","meta");
cd(BasePath)
end

cd(originalPath)