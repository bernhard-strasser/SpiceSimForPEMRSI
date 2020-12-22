%% This script is used to specify data and library file path
% Created by Dr. Chao Ma at MGH
%

[~,machineName]   = system('hostname');
machineName       = machineName(1:end-1);
if strcmp(machineName,'boson')
    dataPath      = '/storage0/home/chao/data/124_pahntom_bs_20170828/';
    libPath       = '/storage0/home/chao/code/support/';
elseif strcmp(machineName,'positrons.mgh.harvard.edu')
    dataPath      = '/home/chao/boson/data/124_pahntom_bs_20170828/';
    libPath       = '/home/chao/boson/code/support/';
elseif strcmp(machineName,'chaomac.local')
    dataPath      = [];
    rawDataPath   = [];
    libPath       = '/Users/chaoma/Dropbox/mgh001273/code/support/';
elseif strcmp(machineName,'positrons.mgh.harvard.edu')

elseif strcmp(machineName,'mgh001273.partners.org')
   
end
% addpath([libPath,'utility/']);
% addpath([libPath,'nonrigid_version23/']);
% addpath([libPath,'siemens_io/']);
% samplingPatternPath = [libPath,'sampling_patterns/'];
% addpath(genpath([libPath,'SPICE/']));