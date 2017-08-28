%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vortex: change point tracking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generates data for Figures 2A-B and 2D-E and writes the data to a file
% that can produce the plots in R.

%% LOAD DATA
clc
clear
close all

% More variables
Vortex_load;
Vortex_variables; % create some variables

%% logicals

do.writeR           = false; % write data to R file for plots?

%% Extra variables

% //REDUCED BAYESIAN OBSERVER

x      	= reshape(stim.samples,btrials,nblocks*nsubj)'; % observed data
H       = stim.changeProb(1); % hazard rate
sN      = stim.std(1); % standard deviation normal distribution
low     = 1;
up      = stim.nPositions(1);

% Reduced Bayesian obsession
[B,lr,d,CPP,MC,r] = redBayes_circ(x,H,sN,low,up);

% Belief updates (human) / belief updates model = lr*d
samples             = reshape(stim.samples,btrials,nblocks*nsubj)'; % observed data
bposition           = reshape(data.bucketPosition,btrials,nblocks*nsubj)'; % bucket position

[data.lr, data.delta] = estimateLR(samples,bposition);

% Reshape for plotting
B           = makeLong(B');
d           = makeLong(d');
MC          = makeLong(MC');
RU          = 1-MC;
CPP         = makeLong(CPP');
r           = makeLong(r');
lr          = makeLong(lr');
data.lr     = makeLong(data.lr');
data.delta  = makeLong(data.delta');

% //SET INDEX
% Exclude trials where estimated learning rate is higher than 99% of group
% estimates (extreme outliers)and where delta is equal 0 (see Nassar 2016)
indx    = [data.lr(data.group==1)<prctile(data.lr(data.group==1),99); data.lr(data.group == 0)<prctile(data.lr(data.group==0),99)]& data.delta > 0;

% Learning rates per subject
for s = 1:nsubj
    idx = data.sub == submat(s) & indx;
    learnrate.subj(s) = mean(data.lr(idx));
end

% Normalise confidence
for s = 1:nsubj
    idx = data.sub == submat(s) & indx;
    data.zConf(idx) = zscore(data.confidence(idx));
end

%% Find change points

% Find points where switch has occured
lag             = 4; % how many samples after change?
patternLength   = lag*2+1;

% Track change points of location
useMean = 0; % use median (0)
[confPatterns, cErr, confChangeFull, confSubz] 	= getChange(data.zConf,stim.change',lag,data.sub,useMean,[],indx);
[lrPatterns, lErr, lrChangeFull, lrSubz]        = getChange(shift(data.lr,1,data.sub,data.block),stim.change',lag,data.sub,useMean,[],indx);
[mcPatterns, mErr, mcChangeFull, mcSubz]        = getChange(MC,stim.change',lag,data.sub,useMean,[],indx);
[cppPatterns, cpErr, cppChangeFull, cppSubz]    = getChange(shift(CPP,1,data.sub,data.block),stim.change',lag,data.sub,useMean,[],indx);

%% Learning rate change point

lrPatterns=[lrPatterns,submat']


%%
% 
 if do.writeR
     fnam = fullfile(datapath, 'RPlotData', 'VortexBet_ChangePointSubjs.txt');
    
%     hdr = {'Cpoint1', 'Cpoint2','Cpoint3','Cpoint4', ...
%         'Cpoint5', 'Cpoint6','Cpoint7','Cpoint8', ...
%         'Cpoint9', 'Subj' };
%     
%     fmt = ('%s\t %s\t %s\t %s\t %s\t  %s\t  %s\t  %s\t  %s\t  %s\t    \n');
%     fid = fopen(fnam, 'w');
%     fprintf(fid, fmt, hdr{:});
%     fclose(fid);
%     
%     data=[ lrPatterns]
    
     dlmwrite(fnam,lrPatterns,'-append','delimiter','\t');
 end

