function [data, stim] = Vortex_initialise(varargin)
% initialise all variables
% function [data, stim] = Vortex_initialise([ppnr], [ntrials], [nblocks])

%% DEFAULT VALUES
optargs = {99 400 10};

% Now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
specif = find(~cellfun(@isempty,varargin)); % find position of specified arguments
[optargs{specif}] = varargin{specif};

% Place optional args in memorable variable names
[ppnr, ntrials, nblocks] = optargs{:};

%% Initialise all necessary variables

% Variables for data
data.sub            = -99*ones(ntrials,1);
data.trial          = [1:ntrials]';
data.block          = zeros(ntrials,1);

data.confidence     = -1;
data.bucketPosition = -1;
data.bet            = -1;
data.hit            = -1;
data.bucketRT       = -1;
data.betRT          = -1;
data.trialReward    = -1;
data.blockReward    = -1;
data.totReward      = -1;

%% Stimulus variables

stim.std         	= -1;
stim.changeProb     = -1;
stim.change         = -1;
stim.steps          = -1;
stim.diameter       = -1;
stim.dotSize        = -1;
stim.bucketRange    = [0 0];
stim.bucketDepth    = -1;
stim.bucketSize     = -1;

stim.nPositions     = 0; % number of possible end points on circle
stim.allCoords      = nan;
stim.xPos           = nan;
stim.yPos           = nan;
stim.refLocation    = nan;
stim.samples0       = nan;
stim.samples        = nan;
stim.betRange       = [0 0];
stim.varChange      = -1;

end