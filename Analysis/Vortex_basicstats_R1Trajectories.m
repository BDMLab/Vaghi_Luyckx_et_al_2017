%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vortex: basic statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Does basic checks for equality between groups.

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
do.saving_plot    	= false; % save plots?
do.writeR           = false; % write data to R file for plots?

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
%
% % //Learning rates
keySet      = submat;
valueSet    = str2num(participant.group);
mapObj      = containers.Map(keySet, valueSet);
% % %

%% Check equal number of change points.
%  Mean per subject
for s = 1:nsubj
    idx = data.sub == submat(s) & indx;
    change.NumberChange_subj(s) = sum(stim.change(idx));
end
% Mean per group
for g = 1:ngroups
    idx = part.group == groups(g);
    change.NumberChange_group(g)        = mean(change.NumberChange_subj(idx));
    change.NumberChange_errgroup(g) 	= std(change.NumberChange_subj(idx))/sqrt(length(change.NumberChange_subj(idx)));
end

[h,p,ci,stats] = ttest2((change.NumberChange_subj(part.group == groups(1))), (change.NumberChange_subj(part.group == groups(2)))) % no stats diff
%

%% Check mean number of trials in between change points
%  Mean per subject

for s = 1:nsubj
    idx = data.sub == submat(s) & indx;
    trials= stim.change(idx);
    trialsChange= find(trials==1); %find ind of trials where change point happens
    diffChange=diff(trialsChange); %compute difference in number of trials in between
    
    change.NumberStable_subj(s)=mean(diffChange);
end

% Mean per group
for g = 1:ngroups
    idx = part.group == groups(g);
    change.group(g)      = mean(change.NumberStable_subj(idx));
    change.err_group(g) 	= std(change.NumberStable_subj(idx))/sqrt(length(change.NumberStable_subj(idx)));
end
[h,p,ci,stats] = ttest2((change.NumberStable_subj(part.group == groups(1))), (change.NumberStable_subj(part.group == groups(2)))) % no stats diff

%% Foucus on position of where samples land

%reshape change point
changes             = reshape(stim.change,btrials,nblocks*nsubj)'; % observed data

diff=0*samples;
rowz    = size(samples,1);
columnz = size(samples,2);

%% Compute diff in sample position from one trial to the next
for b = 1:rowz
    for t = 1:columnz
        if t < columnz
            diff(b,t)    = (diffcirc(samples(b,t),samples(b,t+1)));
        end
    end
end

% Make long both diff2 and changes
diff        = makeLong(diff');
changes     = makeLong(changes');

diff        = [0 ; diff]; %First trial of each subject diff is always 0, add zero on top
% to have diff matching changepoints

%% check for change points
%Exclude trials where change point are not happening
indxC=(changes==1);

%Compute for each subject
%Mean of values for when a change point happens=
%Mean of distance of sample from previous trial when it was stable
for s = 1:nsubj
    idx = data.sub == submat(s) & indx & indxC;
    pos.change(s,1) = mean(diff(idx));
    pos.change(s,2) = submat(s);
    pos.change(s,3) = mapObj(submat(s));
    
end


% % Same per group
for g = 1:ngroups
    idx = part.group == groups(g);
    pos.groupChange(g)      = mean(pos.change(idx));
    pos.grouperrChange(g) 	= std(pos.change(idx))/sqrt(length(pos.change(idx)));
end

[h,p,ci,stats] = ttest2((pos.change(part.group == groups(1))), (pos.change(part.group == groups(2)))) % no stats diff

%% check for intervals where no change points
%Exclude trials where change point are  happening
indxC=(changes==0)

%Compute for each subject
%Mean of values for when a change point happens=
%Mean of distance of sample from previous trial when it was stable
for s = 1:nsubj
    idx = data.sub == submat(s) & indx & indxC;
    pos.Nochange(s,1) = mean(diff(idx));
    pos.Nochange(s,2) = submat(s);
    pos.Nochange(s,3) = mapObj(submat(s));
    
end

% % Same per group
for g = 1:ngroups
    idx = part.group == groups(g);
    pos.groupNoChange(g)      = mean(pos.Nochange(idx));
    pos.grouperrNoChange(g) 	= std(pos.Nochange(idx))/sqrt(length(pos.Nochange(idx)));
end

[h,p,ci,stats] = ttest2((pos.Nochange(part.group == groups(1))), (pos.Nochange(part.group == groups(2)))) % no stats diff

