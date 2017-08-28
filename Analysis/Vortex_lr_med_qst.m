%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vortex_var: basic statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute association between medication and learning rate in patients 
%% LOAD DATA

clc
clear
close all

% More variables
Vortex_load;
Vortex_variables; % create some variables

%% Define more stuff

signif          = 1;
plotting        = 1;
subplotting     = 1;
saving_plot     = 0;
example_subj    = 1;
example_trial   = 1;

%% Logicals

do.plotting      	= true; % plot?
do.signif        	= true; % do significance testing?
do.saving_plot    	= false; % save plots?
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

% //Learning rates

% Learning rates per subject
for s = 1:nsubj
    idx = data.sub == submat(s) & indx;
    learnrate.subj(s) = mean(data.lr(idx));
end

nbins       = 3;
quantSpace  = [min(learnrate.subj) quantile(learnrate.subj,nbins-1) max(learnrate.subj)];
catlr=discretize(learnrate.subj, quantSpace, 'IncludedEdge', 'right');

for g = 1:ngroups
    
    idx = part.group == groups(g);
    [table,chi2,p] = crosstab(catlr(idx),quest.MED(idx))
    
end

for g = 1:ngroups
    
    idx = part.group == groups(g);
    [r,p] = corr(learnrate.subj(idx)',quest.MED(idx), 'type', 'Spearman')
end

for g = 1:ngroups
    
    idx = part.group == groups(g);
    [r,p] = corr(learnrate.subj(idx)',quest.BIS(idx), 'type', 'Pearson')
end
