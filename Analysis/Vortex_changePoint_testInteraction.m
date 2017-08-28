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

% Logicals

do.writeR           = true; % write data to R file for plots?

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
    data.zlr(idx) = zscore(data.lr(idx));
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
[lrPatterns, lErr, lrChangeFull, lrSubz]        = getChange(shift(data.zlr,1,data.sub,data.block),stim.change',lag,data.sub,useMean,[],indx);


zLR_CPP= [lrPatterns, unique(lrSubz), repmat(1,numel(unique(lrSubz)),1)];
zCF_CPP= [confPatterns, unique(confSubz), repmat(2,numel(unique(confSubz)),1)];

zLRCF_CPP= [zLR_CPP; zCF_CPP];

if do.writeR
    filename = fullfile(datapath, 'RPlotData', 'VortexBet_ZVALUES_LRCF_CPP.txt');
    dlmwrite(filename, zLRCF_CPP)
end



%% Learning rate change point

for g = 1:ngroups
    
    idx = part.group == groups(g);
    
    ydat    = nanmean(lrPatterns(idx,:));
    edat    = nanstd(lrPatterns(idx,:))/sqrt(size(lrPatterns(idx,:),1));
    
    HumanLR.group(g,1)=grouplab(g);
    HumanLR.ydat(g,:)=[ydat];
    HumanLR.edat(g,:)=[edat];
    
end

%% Confidence

for g = 1:ngroups
    
    idx = part.group == groups(g);
    
    %ydat    = -1.*nanmean(confPatterns(idx,:)); %this to have reported uncertainty
    ydat    = nanmean(confPatterns(idx,:)); %this to have confidence
    edat    = nanstd(confPatterns(idx,:))/sqrt(size(confPatterns(idx,:),1));
    
    HumanCF.group{g,1}  = cellstr(grouplab(g));
    HumanCF.ydat(g,:)   = [ydat];
    HumanCF.edat(g,:)   = [edat];
end

if do.writeR
    fnam = fullfile(datapath, 'RPlotData', 'VortexBet_ZVALUES_ChangePoint.txt');
    
    hdr = {'OCD_LRydat', 'CTL_LRydat','OCD_LRedat','CTL_LRedat', ...
        'OCD_CFydat', 'CTL_CFydat','OCD_CFedat','CTL_CFedat'};
    
    fmt = ('%s\t %s\t %s\t %s\t %s\t  %s\t  %s\t  %s\t   \n');
    fid = fopen(fnam, 'w');
    fprintf(fid, fmt, hdr{:});
    fclose(fid);
    
    data=[ HumanLR.ydat', HumanLR.edat', ...
        HumanCF.ydat', HumanCF.edat']
    
    dlmwrite(fnam,data,'-append','delimiter','\t');
end
