%%%%%%%%%%%%%%%%%%%%%%
%% Vortex: regressions
%%%%%%%%%%%%%%%%%%%%%%

% Runs the regressions presented in Figure 2C and 2F and writes a datafile
% to produce the plots in R.

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
% estimates (extreme outliers)
indx    = [data.lr(data.group==1)<prctile(data.lr(data.group==1),99); data.lr(data.group == 0)<prctile(data.lr(data.group==0),99)] & data.delta > 0;

% //Confidence

% Normalise confidence
for s = 1:nsubj
    idx = data.sub == submat(s) & indx;
    data.zConf(idx) = zscore(data.confidence(idx));
    data.zBPos(idx) = zscore(data.bucketPosition(idx));
    data.zDiffConf(idx) = [0; abs(diff(-1*zscore(data.confidence(idx))))];
    data.zDiffBPos(idx) = [0; abs(diff(data.bucketPosition(idx)))];
end

%% Possible predictors

X1      = abs(d); % absolute prediction error (model)
X2      = CPP; % change point probability (model)
X3      = RU; % relative uncertainty (model)
X4      = data.hit; % hit or miss?
X5      = data.zDiffBPos'; %absolute difference in bucket update 

shiftX1 = shift(X1,1,data.sub,data.block);
shiftX2 = shift(X2,1,data.sub,data.block);
shiftX4 = shift(X4,1,data.sub,data.block);

%% Regression: confidence

xdat        = [shiftX1 shiftX2 (1-shiftX2).*X3 shiftX4];

ydat        = data.confidence;
rtype       = 'normal';
ptype       = 'linear';
normalise   = 1;

regrlab     = {'|PE|','CPP','(1-CPP)*(1-MC)','hit'};

% Normalise confidence
for s = 1:nsubj
    idx = data.sub == submat(s) & indx;
    ydat(idx) = zscore(ydat(idx));
end

param.subj  = [];
Beta        = [];
G           = [];

for g = 1:ngroups
    
    % Omit first trial (doesn't have subsequent update, so can't be estimated)
    indx1 = ones(ttrials,1);
    indx1(1:btrials:ttrials) = 0;
    
    idx = data.group == groups(g) & indx;
    
    % Regression
    [out, stats, mdl] = regr_fitglm(xdat(idx,:),ydat(idx),data.sub(idx),rtype,ptype,0,normalise,[],indx1(idx));
    
    % Get number of regressors used
    nregr   = size(out.Betas,2)-1;
    
    %Build file for Betas to be plotted
    Beta    = [Beta;out.Betas];
    G       = [G;repmat(groups(g),size(out.Betas,1),1)];
    
end

% Test difference between groups
if do.signif
    
    clear p h
    for r = 1:nregr        
        [p(r),h(r),Cstats(r)] = ranksum(Beta(G==1,r+1),Beta(G==0,r+1));        
    end
    
end

%% Save Betas for Regression confidence

if do.writeR
    fnam    = fullfile(datapath, 'RPlotData', 'VortexBet_RegConfidence.txt');
    hdr     = {'int', '|PE|','CPP','(1-CPP)*(1-MC)','hit', 'Group'};
    Reg     =[Beta, G];
    
    
    fmt = ('%s\t %s\t %s\t %s\t %s\t %s\t \n');
    % fmt(end:end+1) = '\r';
    fid = fopen(fnam, 'w');
    fprintf(fid, fmt, hdr{:});
    fclose(fid);
    dlmwrite(fnam,Reg,'-append','delimiter','\t');
end

%% Regression for learning rate

xdat    = [X1 X1.*X2 X1.*(1-X2).*X3 X1.*X4];
ydatmod = [.5.*X1 lr.*X1]; % fixed and optimal
ydathum = data.lr.*abs(data.delta);
rtype   = 'normal';
ptype   = 'linear';
normalise = 0;

nregr   = size(xdat,2);

regrlab = {'|PE|','CPP','(1-CPP)*(1-MC)','hit'};

% Omit final trial (doesn't have subsequent update, so can't be estimated)
indx1 = ones(ttrials,1);
indx1(1:btrials:ttrials) = 0;

param.subj  = [];
Beta        = [];
G           = [];

for g = 1:ngroups
    
    idx     = data.group == groups(g) & indx;
    
    % Regression human subjects (don't normalise regressors!)
    [out, stats, mdl]  = regr_fitglm(xdat(idx,:),ydathum(idx),data.sub(idx),rtype,ptype,0,normalise,[],indx1(idx));
    
    % Save parameters
    param.subj = [param.subj; out.Betas];
    
    % Get number of regressors used
    nregr   = size(out.Betas,2)-1;
    
    %Build file for Betas to be plotted
    Beta            =[Beta;out.Betas];
    G               =[G;repmat(groups(g),size(out.Betas,1),1)];
end

%% Save Betas for Learning Rate

if do.writeR
    
    fnam    = fullfile(datapath, 'RPlotData', 'VortexBet_RegLearningRate.txt');
    hdr     = {'int', '|PE|','CPP','(1-CPP)*(1-MC)','hit', 'Group'};
    Reg     = [Beta, G];
    
    
    fmt = ('%s\t %s\t %s\t %s\t %s\t %s\t \n');
    % fmt(end:end+1) = '\r';
    fid = fopen(fnam, 'w');
    fprintf(fid, fmt, hdr{:});
    fclose(fid);
    dlmwrite(fnam,Reg,'-append','delimiter','\t');
end

%% Regression: diff conf ~ diff update

xdat        = X5;

ydat        = data.zDiffConf;
rtype       = 'normal';
ptype       = 'linear';
normalise   = 1;

regrlab = {'|diff(Update)|'};

param.subj  = [];
Beta        = [];
G           = [];

for g = 1:ngroups
    
    % Omit first trial (doesn't have subsequent update, so can't be estimated)
    indx1 = ones(ttrials,1);
    indx1(1:btrials:ttrials) = 0;
    
    idx = data.group == groups(g) & indx;
    
    % Regression
    [out, stats, mdl] = regr_fitglm(xdat(idx,:),ydat(idx),data.sub(idx),rtype,ptype,0,normalise,[],indx1(idx));
    
    % Get number of regressors used
    nregr   = size(out.Betas,2)-1;
    
    % Save parameters
    param.subj = [param.subj; out.Betas];
        
    %Save stats
    reg_RegDISC{g} 	= stats;
    out_RegDISC{g} 	= out;
    
    %Build file for Betas to be plotted
    Beta            =[Beta;out.Betas];
    G               = [G;repmat(groups(g),size(out.Betas,1),1)];
    
end

%% Get OCI scores in the correct order as order of the parameters

if do.writeR
    
    idxCTL  = strfind(quest.Group, 'Control');
    idxCTL  = find(not(cellfun('isempty', idxCTL)));
    quest.OCI(idxCTL,2)=1;
    
    idxOCD  = strfind(quest.Group, 'Patient');
    idxOCD  = find(not(cellfun('isempty', idxOCD)));
    quest.OCI(idxOCD,2)=0;
    
    
    OCI_OCD = quest.OCI(idxOCD,:);
    OCI_CTL = quest.OCI(idxCTL,:);
    
    OCI     =[ OCI_OCD; OCI_CTL];
    
    
    %% %% Save Betas for Discrepancy plus OCI scores
    
    fnam    = fullfile(datapath, 'RPlotData', 'VortexBet_RegDiscrepancy.txt');
    hdr     = {'int', 'diff|Update|', 'Group', 'OCI', 'Groupcheck'};
    Reg     = [Beta, G, OCI];
    
    
    fmt = ('%s\t %s\t %s\t %s\t %s\t \n');
    % fmt(end:end+1) = '\r';
    fid = fopen(fnam, 'w');
    fprintf(fid, fmt, hdr{:});
    fclose(fid);
    dlmwrite(fnam,Reg,'-append','delimiter','\t');
end