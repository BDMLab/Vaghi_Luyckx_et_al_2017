%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vortex: basic statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generates the data for plots of 1B and 1C and writes a datafile to produce 
% the plots in R.

%% LOAD DATA
clc
clear
close all

% More variables
Vortex_load;
Vortex_variables; % create some variables

%% Logicals

do.plotting      	= true; % plot?
do.signif        	= true; % do significance testing?
do.saving_plot    	= true; % save plots?
do.writeR           = true; % write data to R file for plots?

%% Set variables

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

[data.lr, data.delta, data.bp] = estimateLR(samples,bposition);

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
data.bp= makeLong(data.bp');
% //SET INDEX
% Exclude trials where estimated learning rate is higher than 99% of group
% estimates (extreme outliers)and where delta is equal 0 (see Nassar 2016)
indx    = [data.lr(data.group==1)<prctile(data.lr(data.group==1),99); data.lr(data.group == 0)<prctile(data.lr(data.group==0),99)]& data.delta > 0;


% //Learning rates
keySet      = submat;
valueSet    = str2num(participant.group);
mapObj      = containers.Map(keySet, valueSet);

% Normalise confidence
for s = 1:nsubj
    idx = data.sub == submat(s) & indx;
    data.zConf(idx) = zscore(data.confidence(idx));
end

%Normalise LR 
for s = 1:nsubj
    idx = data.sub == submat(s) & indx;
    data.zlr(idx) = zscore(data.lr(idx));
end

%% Overal z values learning rate and confidence 
% z Learning rates  and z Confidence per subject
for s = 1:nsubj
    idx = data.sub == submat(s) & indx;
    learnratesubj(s,1) = mean(data.zlr(idx));
    learnratesubj(s,2) = 1
    learnratesubj(s,3) = submat(s);
    learnratesubj(s,4) = mapObj(submat(s)); 
end


% z Learning rates  and z Confidence per subject
for s = 1:nsubj
    idx = data.sub == submat(s) & indx;
    confsubj(s,1) = mean(data.zConf(idx));
    confsubj(s,2) = 2
    confsubj(s,3) = submat(s);
    confsubj(s,4) = mapObj(submat(s)); 
end

testint=[learnratesubj; confsubj]
 
% Store learning rate values into a txt file for R plots

if do.writeR
    filename = fullfile(datapath, 'RPlotData', 'VortexBet_LR_CONF_ZSCORES.txt');
    dlmwrite(filename, testint);
end

%% time windows of 30 trials 
time=(0:30:300);
time(1)=[];

T=[];
learnrate_time=[]; 

for s = 1:nsubj
    
    for t=1:numel(time)
    idx = data.sub == submat(s) & indx & data.trial<=time(t);
    T(t,1) = mean(data.zlr(idx));
    T(t,2) = t;
    T(t,3) = 1
    T(t,4) = submat(s);
    T(t,5) = mapObj(submat(s)); 
    end
    learnrate_time=[learnrate_time;T]
end

T=[]
conf_time=[]

for s = 1:nsubj
    
    for t=1:numel(time)
    idx = data.sub == submat(s) & indx & data.trial<=time(t);
    T(t,1) = mean(data.zConf(idx));
    T(t,2) = t;
    T(t,3) = 2;
    T(t,4) = submat(s);
    T(t,5) = mapObj(submat(s)); 
    end
    conf_time=[conf_time;T]
end

testint_time=[learnrate_time; conf_time]

if do.writeR
    filename = fullfile(datapath, 'RPlotData', 'VortexBet_LR_CONF_ZSCORE_TIME.txt');
    dlmwrite(filename, testint_time);
end

