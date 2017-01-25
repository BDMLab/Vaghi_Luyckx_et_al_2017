function [data, stim, time, aborted] = Vortex_trial(randomise, practice, rotaryConnected, varargin)
% function [data, stim, time, aborted] = Vortex_trial(randomise, practice, rotaryConnected, [ppnr],[ntrials],[nblocks])

%% DEFAULT VALUES

optargs = {99 400 10};

% Now put these defaults into the valuesToUse cell array,
% and overwrite the ones specified in varargin.
specif = find(~cellfun(@isempty,varargin)); % find position of specified arguments
[optargs{specif}] = varargin{specif};

% Place optional args in memorable variable names
[ppnr, ntrials, nblocks] = optargs{:};

%% Setup

% Psychtoolbox defaults
PsychDefaultSetup(2);

% Initialize SoundDriver
InitializePsychSound(1);

% Skip sync tests
Screen('Preference', 'SkipSyncTests', 1);

% Debug
Screen('Preference', 'VisualDebugLevel', 1);

% Put cursor in command window
commandwindow;

%% Keys
KbName('UnifyKeyNames');
space   = KbName('space');
esc  	= KbName('Escape');

alwaysKeys = ([esc, space]);

% Keys with/without rotary
if rotaryConnected == 0
    left            = KbName('LeftArrow');
    right           = KbName('RightArrow');
    noRotaryKeys    = [left, right];
else
    noRotaryKeys = [];
end

RestrictKeysForKbCheck([alwaysKeys, noRotaryKeys]);

%% Data path to save temporary variables

datapath = cd;
datafolder = fullfile(datapath,'/Vortex_data/');

%% Initialise data and setup

[data, stim, time] = Vortex_setup(randomise,ppnr,ntrials,nblocks);

btrials = ntrials/nblocks;

%% Functions

% Make spaced intervals from vector values
makeLinspace = @(x) x(1):x(2);

% Translate points into money
points2money = @(points) round(3 + 4*(points/((1-stim.changeProb)*10*ntrials)));

%%
try
    
    %% Trial variables (based on PTB tutorials by Peter Scarfe)
    
    % //SCREEN variables
    
    screens         = Screen('Screens');
    screenNumber    = max(screens);
    
    % COLOUR variables
    
    % Basic colours
    col.white           = WhiteIndex(screenNumber); % Define black and white (white will be 1 and black 0). This is because
    col.gradient        = repmat(col.white*.9,1,3);
    col.black           = BlackIndex(screenNumber); % luminace values are (in general) defined between 0 and 1.
    
    % Background colours
    col.background      = col.white / 2;
    col.defaultFrame    = col.white .* .9;
    col.green           = [0 1 0] .*col.white;
    col.red             = [1 0 0] .*col.white;
    
    % Bucket colours
    bucketcolz          = [zeros(stim.bucketRange(1),3); autumn(stim.bucketRange(2)-stim.bucketRange(1))];
    col.bucket          = bucketcolz;
    
    % Money colours
    col.money           = [1 .9 0]*col.white;
    
    % Coin icon colour
    col.iconFront  	= [.99; 0.89; 0]*col.white;
    col.iconBack  	= [.98; 0.65; 0]*col.white;
    
    % Open window and window data
    [w, windowRect]                 = PsychImaging('OpenWindow', screenNumber, col.background); % Open an on screen window and color it grey
    [screenXpixels, screenYpixels]  = Screen('WindowSize', w);     % Get the size of the on screen window in pixels
    [scr.xCenter, scr.yCenter]      = RectCenter(windowRect); % Get the centre coordinate of the window in pixels
    
    Screen('BlendFunction', w, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    
    % //TEXT VARIABLES
    
    leftMargin      = 50;
    rightMargin     = leftMargin;
    topMargin       = 50;
    
    titleSize       = 60;
    textSize        = 26;
    scoreSize       = 30;
    standardFont    = 'Calibri';
    scoreFont       = 'Century Gothic';
    
    Screen('TextFont', w, standardFont); % Font
    
    % //STIMULUS variables
    
    % Fixation dot
    scr.fixRad = 9; % radius of the fixation circle
    
    % Frame variables
    scr.penWidthPixels 	= 10; % Pen width for the frame
    diameter            = stim.diameter; % diameter of big circle
    
    stimRect            = [0 0 diameter diameter]; % Frame for stimulus presentation
    circleXpos          = screenXpixels * .5;
    circleYpos          = screenYpixels * .5;
    scr.circleCoord  	= CenterRectOnPointd(stimRect, circleXpos, circleYpos);
    
    % Bucket variables
    bucketSpeed         = 2;
    bDiameter           = diameter + stim.bucketDepth;
    
    bucketRect          = [0 0 bDiameter bDiameter]; % Frame for stimulus presentation
    bucketXpos          = screenXpixels * .5;
    bucketYpos          = screenYpixels * .5;
    scr.bucketCoord  	= CenterRectOnPointd(bucketRect, bucketXpos, bucketYpos);
    
    % Money variables
    scr.moneyRad        = stim.dotSize;
    
    Screen('TextFont', w, scoreFont);
    Screen('TextSize', w,scoreSize);
    line1               = sprintf('%010s',num2str(0));
    scoreBounds         = Screen('TextBounds', w, line1, 0, 0);
    
    scr.boxPenWidth     = 2;
    boxMargin           = 10;
    boxHeight           = scoreSize + boxMargin;
    boxWidth            = scoreBounds(3) + boxMargin*2;
    scoreRect           = [0 0 boxWidth boxHeight]; % Frame for stimulus presentation
    scr.scoreXpos      	= screenXpixels - boxWidth/2 - rightMargin;
    scr.scoreYpos      	= textSize + 40;
    scr.moneyBoxCoord  	= CenterRectOnPointd(scoreRect, scr.scoreXpos, scr.scoreYpos); % frame box
    
    scr.coinIconRad     = 9;
    coinOffset          = 3;
    scr.coinXpos        = screenXpixels - boxWidth - rightMargin - 30;
    scr.coinIconCoord(:,1)  = [scr.coinXpos + coinOffset; scr.scoreYpos + coinOffset/2];
    scr.coinIconCoord(:,2)  = [scr.coinXpos; scr.scoreYpos];
    
    % Confidence bar
    confBarHeight       = 20;
    confBarWidth        = 500;
    confBarRect         = [0 0 confBarWidth confBarHeight];
    scr.confBarXpos     = screenXpixels * .5;
    scr.confBarYpos     = screenYpixels * .9;
    scr.confBarCoord    = CenterRectOnPointd(confBarRect, scr.confBarXpos, scr.confBarYpos);
        
    confIndHeight       = confBarHeight;
    confIndWidth        = confBarWidth/100;
    scr.confIndRect     = [0 0 confIndWidth confIndHeight];   
    
    % Translate position to confidence
    pos2conf    = @(pos) round(round(((pos - scr.confBarCoord(1))/(scr.confBarCoord(3) - scr.confBarCoord(1))).*100).*99/100 +1);
    
    % //AUDIO variables
    
    nrchannels = 2;
    freq = 48000;    
    repetitions = 1;    
    startCue = 0;
    waitForDeviceStart = 1;
    pahandle = PsychPortAudio('Open', [], 1, 0, freq, nrchannels);    
    PsychPortAudio('Volume', pahandle, 0.5);
    inc1 	= MakeBeep(453*3/4, time.beep, freq); 
    inc2 	= MakeBeep(453/2, time.beep, freq);
    cor1    = MakeBeep(800/2, time.beep, freq);
    cor2    = MakeBeep(800, time.beep, freq);
     
    %% Instructions before experiment
    
    % Hide the mouse cursor
    HideCursor;
    
    % Start of real experiment
    time.expStart = GetSecs;
    
    if practice == 1
        % Instructions before practice
        instructions    = 'beforepractice'; Vortex_instructions;
    elseif practice == 0
        % Instructions after practice
        instructions    = 'afterpractice'; Vortex_instructions;
    end
    
    %% Actual trials
    
    for t = 1:ntrials
        
        %% Start of a block
        
        % Instructions before each new block
        if mod(t,btrials) == 1
            instructions = 'startblock'; Vortex_instructions;
        end
        
        %% Start of trial
        
        if practice == 0
            disp(num2str(t));   % display trial number
        end
        
        % Print big circle
        Screen('TextSize', w,scoreSize);
        Screen('TextStyle',w,0);
        Screen('TextFont', w, scoreFont);
        SetMouse(scr.xCenter, scr.yCenter); % put the cursor in the middle of the screen
        time.trialStart = GetSecs;
        
        %% Bucket positioning
        move        = 0; % shift of bucket
        confirm     = 0; % confirm final position before shooting dot
        %selapsed    = 0;
        
        while confirm == 0 %&& selapsed < time.deadline
            
            moved       = 0; % has bucket moved?
            
            % Stuff to draw
            DefaultFrames(w,scr,col,data.bucketPosition(t,1),stim.bucketSize(t,1),data.hit(t,1));
            MoneyBox(w,scr,col,data.blockReward(t,1));
            Screen('Flip',w);
            
            [kdown, ~, codes] = KbCheck;
            if kdown
                if codes(esc)
                    
                    % Save data
                    if practice == 0
                        tmpname = sprintf('datatmp_%d.mat',ppnr);
                        save(fullfile(datafolder,tmpname),'data','stim','time');
                        RestrictKeysForKbCheck([]);
                        ShowCursor;
                    end
                    
                    aborted = 1;
                    % Close the audio device
                    PsychPortAudio('Close', pahandle);
                    Screen('CloseAll');
                    ShowCursor;
                    return;
                    
                elseif codes(space)
                    confirm = 1;
                    data.bucketRT(t,1) = GetSecs - time.trialStart;
                else
                    
                    % Move bucket
                    if rotaryConnected == 0
                        if codes(left)
                            move    = bucketSpeed;
                            moved   = 1;
                        elseif codes(right)
                            move    = -1*bucketSpeed;
                            moved   = 1;
                        end
                    end
                end
            end
            
            % Mouse recording
            if rotaryConnected == 1
                [mouseX,~,buttons] = GetMouse;
                if any(buttons) == 1
                    confirm = 1;
                end
                move = 0;
                if mouseX ~= scr.xCenter
                    move = bucketSpeed*round((mouseX - scr.xCenter)/10);
                    SetMouse(scr.xCenter, scr.yCenter); % put the cursor in the middle of the screen
                    moved = 1;
                end
            end
            
            % Position of bucket
            if moved == 1
                data.bucketPosition(t,1) = data.bucketPosition(t,1) + move;
            end
        end
        
        % Wait until key is released
        if KbCheck
            KbWait([], 1); 
        end
        
        %% Betting on hit
        
        WaitSecs(time.ISI);
                
        betStart = GetSecs;
        
        finished    = 0;
        selapsed    = 0;
        newX        = (scr.confBarCoord(3) - scr.confBarCoord(1)) * data.bet(t,1)/100 + scr.confBarCoord(1);
        prevX       = newX;
        
        while finished == 0

            moved = 0; % has bucket moved?
            
            [kdown, ~, codes] = KbCheck;
            if kdown
                if codes(esc)
                    
                    % Save data
                    if practice == 0
                        tmpname = sprintf('datatmp_%d.mat',ppnr);
                        save(fullfile(datafolder,tmpname),'data','stim','time');
                        RestrictKeysForKbCheck([]);
                        ShowCursor;
                    end
                    
                    aborted = 1;
                    % Close the audio device
                    PsychPortAudio('Close', pahandle);
                    Screen('CloseAll');
                    ShowCursor;
                    return;
                    
                elseif codes(space) && selapsed > time.RSI
                    finished = 1;
                    data.betRT(t,1) = GetSecs - betStart;
                else
                    
                    % Move numbers
                    if rotaryConnected == 0
                        if codes(left)
                            move    = -5;
                            moved   = 1;
                        elseif codes(right)
                            move    = +5;
                            moved   = 1;
                        end
                    end
                end
            end
            
            % Mouse recording
            if rotaryConnected == 1
                [mouseX] = GetMouse;
                move = 0;
                if mouseX ~= scr.xCenter
                    move = round((mouseX - scr.xCenter));
                    SetMouse(scr.xCenter, scr.yCenter); % put the cursor in the middle of the screen
                    moved = 1;
                end
            end
            
            % Calculate new position of pointer
            if moved == 1
                newX = prevX + move;                
                if newX < scr.confBarCoord(1) || newX > scr.confBarCoord(3)
                    newX = prevX;
                end
                prevX = newX;
            end

            % Visuals
            DefaultFrames(w,scr,col,data.bucketPosition(t,1),stim.bucketSize(t,1),data.hit(t,1));
            MoneyBox(w,scr,col,data.blockReward(t,1));
            ConfidenceIndicator(w,scr,col,newX);
            
            % Flip
            Screen('Flip',w);   
            
            % Check time elapsed
            selapsed    = GetSecs - betStart;
        end
                
        %% Warning for shoot
        
        DefaultFrames(w,scr,col,data.bucketPosition(t,1),stim.bucketSize(t,1),data.hit(t,1));
        dotCoords = [stim.xTrajectory(stim.samples(t),1), stim.yTrajectory(stim.samples(t),1)];
        MoneyDots(w,scr,col,dotCoords);
        MoneyBox(w,scr,col,data.blockReward(t,1));
        %BetIndicator(w,scr,col.white,data.bet(t,1));
        
        Screen('Flip',w);
        WaitSecs(time.ISI);
        
        %% Calculate if dot hits bucket
        
        % Potential range of coordinates where dot might fall
        alpha           = mod(data.bucketPosition(t,1)-stim.bucketSize(t,1)/2,360)+1; % 1 - 360
        beta            = mod(data.bucketPosition(t,1)+stim.bucketSize(t,1)/2,360)+1; % 1 - 360
        bucketCoords    = coordBucket([stim.xPos stim.yPos],alpha,beta);
        
        % Anticipate location of dot arriving
        if ~isempty(intersect(bucketCoords,[stim.xTrajectory(stim.samples(t),end),stim.yTrajectory(stim.samples(t),end)],'rows'))
            data.hit(t,1)           = 1;
        else
            data.hit(t,1)           = 0;
        end     
        
        data.bet(t,1)           = pos2conf(prevX);      
        data.trialReward(t,1)   = 10*(2*data.hit(t,1)-1);
        
        %% Stimulus presentation
        
        for i = 1:time.dotSpeed
            
            % Flash when dot hits bucket
            if i == time.dotSpeed
                if data.hit(t) == 1
                    col.bucket = repmat(col.white,size(bucketcolz,1),3);
                end
            else
                col.bucket   = bucketcolz;
            end
            
            DefaultFrames(w,scr,col,data.bucketPosition(t,1),stim.bucketSize(t,1),-1);
            dotCoords = [stim.xTrajectory(stim.samples(t),i), stim.yTrajectory(stim.samples(t),i)];
            MoneyDots(w,scr,col,dotCoords);
            MoneyBox(w,scr,col,data.blockReward(t,1));
            %BetIndicator(w,scr,numcol,data.bet(t,1));
            Screen('Flip', w);
            
        end
        
        %% End of trial
        
        col.bucket              = bucketcolz;
        data.blockReward(t,1)   = sum(data.trialReward(data.block == data.block(t),:));
        
        DefaultFrames(w,scr,col,data.bucketPosition(t,1),stim.bucketSize(t,1),data.hit(t,1));
        MoneyBox(w,scr,col,data.blockReward(t,1));
        %BetIndicator(w,scr,numcol,data.bet(t,1));
        
        [time.trialEnd]     = Screen(w,'Flip');
        time.trialDur(t)    = time.trialEnd - time.trialStart;
        
        % Audio feedback
        if data.hit(t) == 1
            PsychPortAudio('FillBuffer', pahandle, [cor1; cor2]);
            PsychPortAudio('Start', pahandle, repetitions, startCue, waitForDeviceStart);
        else
            PsychPortAudio('FillBuffer', pahandle, [inc1; inc2]);
            PsychPortAudio('Start', pahandle, repetitions, startCue, waitForDeviceStart);            
        end

        WaitSecs(time.ITI);
        PsychPortAudio('Stop', pahandle);
        
        % Things to save
        data.bucketPosition(t)      = mod(data.bucketPosition(t),360);
        data.confidence(t,1)        = pos2conf(prevX);
        data.totReward(t,1)         = sum(data.trialReward);
        
        % Maintain bucket position and size for next trial
        if mod(t,btrials) ~= 0
            data.bucketPosition(t+1)    = data.bucketPosition(t);
            data.blockReward(t+1)       = data.blockReward(t);
        end
        
        %% End of block
        
        % Break after have a block with practice
        if mod(t,btrials) == btrials && practice == 1
            WaitSecs(time.EBI);
            % Text on screen
            instructions = 'endblock'; Vortex_instructions;
            break;
        end
        
        if mod(t,btrials) == 0
            
            WaitSecs(time.EBI);
            
            % Text on screen
            instructions = 'endblock'; Vortex_instructions;
            
            if practice == 0
                % Save data from block
                tmpname = sprintf('datatmp_%d.mat',ppnr);
                save(fullfile(datafolder,tmpname),'data','stim','time');
            end
        end
        
    end
    
    %% End of experiment
    
    if practice == 0
        % Translate points into money
        money = points2money(data.totReward(t,1))
        
        % Instructions at end of experiment
        instructions = 'endexp'; Vortex_instructions;
        
        % Get length of experiment
        time.expEnd     = GetSecs;
        time.expDur     = time.expEnd - time.expStart;
        
        % Close all
        ShowCursor;
    end
    
    aborted = 0;
      
    % Close the audio device
    PsychPortAudio('Close', pahandle);
    
    RestrictKeysForKbCheck([]);
    Screen('CloseAll');
    
catch
    
    % Get length of experiment.
    if time.expStart > 0
        time.expEnd     = GetSecs;
        time.expDur     = time.expEnd - time.expStart;
    end
    
    aborted = 1;
       
    % Close the audio device
    PsychPortAudio('Close', pahandle);
    
    Screen('CloseAll');
    RestrictKeysForKbCheck([]);
    ShowCursor;
    
    warning('Something went wrong in the experiment.');
    
    if practice == 0
        tmpname = sprintf('datatmp_%d.mat',ppnr);
        save(fullfile(datafolder,tmpname),'data','stim','time');
    end
    
    %FlushEvents;
    rethrow(lasterror);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% // FUNCTIONS

%% Draw default frames
function DefaultFrames(w,scr,col,center,bucketSize,hit)

% Big circle
Screen('FrameOval', w, col.defaultFrame, scr.circleCoord, scr.penWidthPixels);

% Bucket
leftSide  = center - bucketSize/2;
rightSide = bucketSize;

Screen('FillArc', w, col.bucket(bucketSize,:), scr.bucketCoord, leftSide, rightSide);
Screen('FillArc', w, col.background, scr.circleCoord, leftSide, rightSide);

% Centre dot
if hit == 1
    dotCol = col.green;
elseif hit == 0
    dotCol = col.red;
else
    dotCol = col.defaultFrame;
end

Screen('DrawDots', w, [scr.xCenter, scr.yCenter] , scr.fixRad, dotCol, [], 2);

end

%% Draw money dots
function MoneyDots(w,scr,col,whichDotCoords)
Screen('DrawDots', w, whichDotCoords, scr.moneyRad, col.money, [scr.xCenter, scr.yCenter], 2);
end

%% Get coordinates of all locations in bucket range
function [coords] = coordBucket(allCoords,alpha,beta)

alpha   = floor(alpha);
beta    = ceil(beta);

if alpha > 360
    alpha = 360;
end
if beta > 360
    beta = 360;
end

coords = allCoords(alpha:beta,:);

if isempty(coords)
    ind     = [alpha:360 1:beta];
    coords  = allCoords(ind,:);
end

end

%% Confidence indicator

function ConfidenceIndicator(w,scr,col,newX)

scr.confIndCoord    = CenterRectOnPointd(scr.confIndRect, newX, scr.confBarYpos);
    
Screen('FillRect',w,col.defaultFrame,scr.confBarCoord);
Screen('FillRect',w,col.black,scr.confIndCoord);

line1 = '1';
line2 = '100';
rct1 = CenterRectOnPoint(Screen('TextBounds',w,line1),scr.confBarCoord(1),scr.confBarYpos + 30);
rct2 = CenterRectOnPoint(Screen('TextBounds',w,line2),scr.confBarCoord(3),scr.confBarYpos + 30);
Screen('DrawText',w,line1,rct1(1),rct1(2),col.white);
Screen('DrawText',w,line2,rct2(1),rct2(2),col.white);

end

%% Money indicator

function MoneyBox(w,scr,col,points)

% Draw score
line1 = sprintf('%d',points);
rct = CenterRectOnPoint(Screen('TextBounds',w,line1),scr.scoreXpos,scr.scoreYpos);
Screen('DrawText',w,line1,rct(1),rct(2),col.white);

% Draw coin icon
Screen('DrawDots', w, scr.coinIconCoord, scr.coinIconRad, [col.iconBack col.iconFront], [], 2);

end
