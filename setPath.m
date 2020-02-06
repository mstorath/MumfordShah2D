function setPath
javaaddpath(fullfile(fileparts(which(mfilename)), 'Java', 'bin'))
addpath(fullfile(fileparts(which(mfilename)), '.'))
addpath(fullfile(fileparts(which(mfilename)), '.', 'Auxiliary'))
end

