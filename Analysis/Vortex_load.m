
% Fill in location of the folder. Generates path names for all necessary
% folders. Loads behavioural and questionnaire data.

%% LOAD DATA

%Set path
mainfolder      = '~/Dropbox/Vortex_bet/GitHub'; % fill in location of folder
mypath          = fullfile(mainfolder,'Analysis');
datapath        = fullfile(mainfolder,'Data');
figpath         = fullfile(mainfolder,'Plots');
functionpath    = fullfile(mainfolder,'Functions');

addpath(figpath);
addpath(genpath(datapath));
addpath(genpath(functionpath));
cd(mypath);

% Load experiment data
load(fullfile(datapath,'Behavioural' , 'Vortex_ocd_fulldata.mat'));

% Load questionnaire data
quest = table2struct(readtable(fullfile(datapath,'Questionnaires' , 'VortexBet_DepOCDMeasures.txt')),'ToScalar',true);

quest.YBOCS_TOT = str2double(quest.YBOCS_TOT);
quest.YBOCS_OBS = str2double(quest.YBOCS_OBS);
quest.YBOCS_COM = str2double(quest.YBOCS_COM);
quest.COINS = str2double(quest.COINS);
quest.AES = str2double(quest.AES);
quest.BISBAS = str2double(quest.BISBAS);
quest.MED = str2double(quest.MED);



