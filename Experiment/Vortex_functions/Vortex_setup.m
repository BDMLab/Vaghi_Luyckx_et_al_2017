function [data, stim, time] = Vortex_setup(randomise, varargin)
% [data, stim] = Vortex_setup(randomise, [ppnr],[ntrials], [nblocks])

%% DEFAULT VALUES

optargs = {99 400 10};

% Now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
specif = find(~cellfun(@isempty,varargin)); % find position of specified arguments
[optargs{specif}] = varargin{specif};

% Place optional args in memorable variable names
[ppnr, ntrials, nblocks] = optargs{:};

%% Initialise variables

[data, stim] = Vortex_initialise(ppnr,ntrials,nblocks);

btrials = ntrials/nblocks;

makeLong = @(x) x(:); % function to squelch data

%% Get screen variables

nScreens        = Screen('Screens');
screenNumber    = max(nScreens);
screen          = Screen('Resolution',screenNumber);

%% Time variables

time.framedur       = 1/screen.hz;                  % approximate frame duration
if time.framedur == Inf
    time.framedur = .0167;
    warning('No frame duration registered.');
end

time.deadline       = 2.5;                          % deadline for response
time.ISI            = .15;                          % time between confirm and dot shoot
time.RSI            = .5;                           % time between two key presses
time.ITI            = .8;                           % intertrialinterval
time.EBI            = 1;                            % end block interval
time.dotSpeed       = floor(.5/time.framedur);      % how many frames to reach other side? (seconds -> frames)
time.beep           = time.ITI/2;                   % length of sound feedback

time.expStart       = 0;                    % start of experiment
time.expEnd         = 0;                    % end of experiment
time.expDur         = 0;                    % duration of experiment
time.trialStart     = 0;                    % start of trial
time.trialEnd       = 0;                    % end of trial
time.trialDur       = zeros(ntrials,1);     % duration of trial
time.presEnd        = 0;                    % end of sample presentation

%% Stimulus variables

% Create stimulus data
SD                  = [12];
stim.std            = repmat(repmat(SD,btrials/length(SD),1),nblocks,1); % category of variance

stim.changeProb     = .125;
stim.change         = zeros(ntrials,1);
stim.varChange      = zeros(ntrials,1);

stim.diameter       = 500; % diameter of big circle (every degree = 4 px)
stim.dotSize        = 8; % radius of money dots

minSize             = 1*2;
maxSize             = stim.std(1)*3*2;
stim.bucketRange    = [minSize maxSize]; % min and max size of bucket (in degrees) - 1
stim.bucketDepth    = 50;
stim.bucketSize     = repmat(SD*3,ntrials,1);

stim.steps          = 4; % how much drift (degrees)
stim.nPositions     = 360; % number of possible end points on circle
stim.allCoords      = circle(stim.nPositions,stim.diameter/2 - stim.dotSize/2);
stim.xPos           = stim.allCoords(:,1);
stim.yPos           = -1*stim.allCoords(:,2); % get in direction of screen y axis
for i = 1:stim.nPositions
    stim.xTrajectory(i,:) = linspace(0,stim.xPos(i),time.dotSpeed);
    stim.yTrajectory(i,:) = linspace(0,stim.yPos(i),time.dotSpeed);
end
stim.refLocation    = repmat(180,ntrials,1);

stim.betRange       = [1 100]; % min and max value of coins

%% Data variables

% Variables for data
data.sub            = ppnr*ones(ntrials,1);
data.block          = makeLong(repmat(1:nblocks,btrials,1));
data.bet            = 50*ones(ntrials,1);
data.hit            = -1*ones(ntrials,1);
data.bucketPosition = 90*ones(ntrials,1);
data.confidence     = data.hit;
data.bucketRT       = zeros(ntrials,1);
data.betRT          = data.bucketRT;
data.trialReward    = zeros(ntrials,1);
data.blockReward    = zeros(ntrials,1);

%% Randomise everything

if randomise == 1
    
    %% Change of reference location
    for t = 1:ntrials
        if rand(1) <= stim.changeProb || mod(t,btrials) == 1
            currentLocation = 1 + (stim.nPositions-1).*rand(1);
            stim.change(t)  = 1;
            stim.refLocation(t) = round(currentLocation);
        else
            stim.refLocation(t) = round(stim.refLocation(t-1));
        end      
    end
    
    % Randomise position of bet indicator
    data.bet = randi(50,ntrials,1)+25;
    
end

%% Get sample timing and trajectories

% Draw samples around zero
stim.samples0       = round(bsxfun(@times,stim.std,randn(ntrials,1)));

% Put samples around reference location
stim.samples        = stim.samples0 + stim.refLocation;
stim.samples        = mod(stim.samples,360);
stim.samples(stim.samples == 0) = 360;

end
