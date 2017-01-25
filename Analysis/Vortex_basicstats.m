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
keySet      = submat;
valueSet    = str2num(participant.group);
mapObj      = containers.Map(keySet, valueSet);

% Learning rates per subject
for s = 1:nsubj
    idx = data.sub == submat(s) & indx;
    learnratesubj(s,1) = mean(data.lr(idx));
    learnratesubj(s,2) = submat(s);
    learnratesubj(s,3) = mapObj(submat(s));
    
end

%% Store learning rate values into a txt file for R plots

if do.writeR
    filename = fullfile(datapath, 'RPlotData', 'VortexBet_LR.txt');
    dlmwrite(filename, learnratesubj);
end

%% Additional measures relating to accuracy and confidence and change points

% Hits and points per subject, per block
for s = 1:nsubj
    for b = 1:nblocks
        idx     = data.sub == submat(s) & data.block == b & indx;
        hit.subj_blocks(s,b)    = mean(data.hit(idx));
        point.subj_blocks(s,b)  = mean(data.trialReward(idx));
    end
end

for g = 1:ngroups
    for b = 1:nblocks
        hit.group_blocks(g,b)   = mean(hit.subj_blocks(part.group == groups(g),b));
        point.group_blocks(g,b) = mean(point.subj_blocks(part.group == groups(g),b));
    end
end

% Hits and points per subject
hit.subj        = mean(hit.subj_blocks,2);
hit.group       = mean(hit.group_blocks,2);
hit.err_subj   	= std(hit.subj)/sqrt(nsubj);
point.subj      = mean(point.subj_blocks,2);
point.group     = mean(point.group_blocks,2);
point.err       = std(point.subj)/sqrt(nsubj);

for g = 1:ngroups
    subzH = hit.subj(part.group == groups(g));
    subzP = point.subj(part.group == groups(g));
    hit.err_group(g)   = std(subzH)/sqrt(length(subzH));
    point.err_group(g) = std(subzP)/sqrt(length(subzP));
end

% //Confidence

% Normalise confidence
for s = 1:nsubj
    idx = data.sub == submat(s) & indx;
    data.zConf(idx) = zscore(data.confidence(idx));
end


% Check equal change point per group

for s = 1:nsubj
    idx = data.sub == submat(s) & indx;
    change.subj(s) = sum(stim.change(idx));
end

% Change point per group
for g = 1:ngroups
    idx = part.group == groups(g);
    change.group(g)      = mean(change.subj(idx));
    change.err_group(g) 	= std(change.subj(idx))/sqrt(length(change.subj(idx)));
end
%% Model behaviour

if do.plotting
    
    s   = 1; % subject 1
    b   = 1; % block 1
    pos = [1,2];
    fig = figure; hold on;
    
    % Trial example    
    idx     = data.sub == submat(s) & data.block == b;
    subplot(2,1,pos(b,1)); hold on;
    modplot     = plot(1:btrials,B(idx),'--','Color','k','LineWidth',lnwid);
    sampplot 	= plot(1:btrials,stim.samples(idx),'yo','MarkerFaceColor',[1 0.6 0],'MarkerEdgeColor',[0.8 0.5 0],'MarkerSize',mksz,'LineWidth',1.2);
    xlim([0 75])
    ylim([0 370])
    xlabel('Trial','FontSize',axfntsz);
    ylabel('Degrees','FontSize',axfntsz);
    title('Beliefs of model' ,'FontSize',titlefntsz);
    
    % Model example
    subplot(2,1,pos(b,2)); hold on;
    cppplot     = plot(1:btrials,CPP(idx),'Color',[0.5 0.5 0.5],'LineWidth',lnwid);
    mcplot      = plot(1:btrials,MC(idx),'Color','k','LineWidth',lnwid);
    xlabel('Trial','FontSize',axfntsz);
    ylabel(' a.u.','FontSize',axfntsz);
    title({'Change point probability (CPP) and model confidence (MC)' },'FontSize',titlefntsz);
    hL = legend([cppplot,mcplot],{'CPP','MC'},'FontSize',axfntsz,'LineWidth',lnwid, 'Location', 'best');
    legend boxoff
    xlim([0 75])
    
    if do.saving_plot
        figname=fullfile(mainfolder, 'Plots', 'Plot_ModelBehaviour')
        saveas(fig,figname, 'epsc');
    end
    
end

%% Ttest learning rate, confidence, change point, hit/missed

if do.signif
    
     [h,p,ci,stats] = ttest2((confidence.subj(part.group == groups(1))), (confidence.subj(part.group == groups(2)))) % no stats diff
    
    [h,p,ci,stats] = ttest2((change.subj(part.group == groups(1))), (change.subj(part.group == groups(2)))) % no stats diff
    
    [h,p,ci,stats] = ttest2((hit.subj(part.group == groups(1))), (hit.subj(part.group == groups(2)))) % no stats diff
end


