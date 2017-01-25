%%%%%%%%%%%%%%%%%%%%%%
%% Vortex_main: runner
%%%%%%%%%%%%%%%%%%%%%%

execute = questdlg('Clear all and start experiment?');
%execute = 'Yes';

if strmatch(execute, 'Yes')
    
    clc
    clear all
    close all
    
    %% Set path
    datapath = cd;
    datafolder = fullfile(datapath,'/Vortex_data/');
    
    addpath(genpath(fullfile(datapath,'Vortex_functions/')));
    addpath(datafolder);
    
    %% Create participant structure
    
    argindlg = inputdlg({'Participant number   ','Gender (M/F)','Age','Hand (L/R)','Control'},'',1,{'000','','','R',''});
    if isempty(argindlg)
        return;
    else
        participant = struct;
        participant.name        = upper(argindlg{1});
        participant.gender      = argindlg{2};
        participant.age         = argindlg{3};
        participant.handedness  = argindlg{4};
        participant.group       = argindlg{5};
    end
    
    %% Get OS
    os = computer;
    
    %% Rotary switch connected?
    rotaryConnected = questdlg('Rotary switch connected?');
    switch rotaryConnected
        case 'Yes'
            rotaryConnected = 1;
%             if all(os == 'PCWIN64') % open powermate in windows
%                 !C:\Program Files (x86)\Griffin Technology\PowerMate\PowerMate.exe &
%             end
        case 'No'
            rotaryConnected = 0;
        otherwise
            return;
    end
    
    %% Set keys
    KbName('UnifyKeyNames');
    
    %% What I need to know
    ppnr            = str2num(participant.name);    % participant number
    nblocks         = 4;                            % number of blocks
    ntrials         = 75*nblocks;                	% number of trials
    aborted         = 0;                            % aborted experiment?
    
    randomise   = 1;                            % variable to activate randomisation
    if randomise == 1
        % Seed random number generator
        rng('shuffle');
    end
    
    %% Demo
    practice = 3;
    disp('Demo initiated.')
    [~, ~, ~, aborted] = Vortex_trial(randomise,practice,rotaryConnected,ppnr,4,1);
    
    %% Practice trials
    if aborted == 0
        practice = 1;
        disp('Practice initiated.')
        [~, ~, ~, aborted] = Vortex_trial(randomise,practice,rotaryConnected,ppnr,20,1);
    end
    
    %% Run experiment
    if aborted == 0
        practice = 0;
        disp('Experiment initiated.')
        [data, stim, time, aborted] = Vortex_trial(randomise, practice, rotaryConnected, ppnr, ntrials, nblocks);
    end
    
    if aborted == 0
        %% Save data
        participant.filename = sprintf('Vortex_ppt_%s_%3s.mat',participant.name,datestr(now,'yyyymmddHHMMSS'));
        datafile = fullfile(datapath,'/Vortex_data/',participant.filename);
        save(datafile,'participant','data','stim','time');
        
        %% Send e-mail with data
      %  Vortex_send_email;
    end
    
end

%% EXPLANATION DIFFERENT FUNCTIONS
%
% Vortex_initialise = define all variables for data and stim(ulus) structures,
% very handy to check when something went wrong in the randomisation or
% trial files.
%
% Vortex_setup = give all variables the real values and option to randomise
% everything. Use this to check whether randomisation does what it's
% supposed to do and check underlying distributions.
%
% Vortex_trial = run the real experiment without all the fuzz around it
%
% Vortex_main = file that executes the whole experiment. Here you define the
% actual values of the number of trials etc. in your experiment.













