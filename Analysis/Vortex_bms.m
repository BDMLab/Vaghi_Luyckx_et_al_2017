%%%%%%%%%%%%%%
%% Vortex: BIC
%%%%%%%%%%%%%%

% Runs 4 possible models of the two main regression models presented in the
% Supplementary Material (Table S3 and S4), calculates the BIC for these 
% models and produces the plots S1 and S2.
% Requires function from SPM toolbox (spm_BMS).

%% LOAD DATA

clc
clear
close all

Vortex_load;

% More variables
Vortex_variables; % create some variables

% Add SPM12 toolbox
addpath('~/Documents/MATLAB/spm12')

%% Logicals

do.plotting     = true;
do.saving_plot  = false;

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

%% Possible predictors

X1      = abs(d); % absolute prediction error (model)
X2      = CPP; % change point probability (model)
X3      = RU; % relative uncertainty (model)
X4      = data.hit; % hit or miss?

shiftX1 = shift(X1,1,data.sub,data.block);
shiftX2 = shift(X2,1,data.sub,data.block);
shiftX4 = shift(X4,1,data.sub,data.block);

%% Run different regression models: confidence

modelz  = {[shiftX1],[shiftX1 shiftX2],[shiftX1 shiftX2 (1-shiftX2).*X3],[shiftX1 shiftX2 (1-shiftX2).*X3 shiftX4]};
avBIC   = zeros(size(modelz,2),ngroups);

ydat        = data.confidence;
rtype       = 'normal';
ptype       = 'linear';
normalise   = 1;

regrlab = {'|PE|','CPP','(1-CPP)*(1-MC)','hit'};

% Normalise confidence
for s = 1:nsubj
    idx = data.sub == submat(s) & indx;
    ydat(idx) = zscore(ydat(idx));
end

for m = 1:size(modelz,2)
    
    disp(['Running regression model ' num2str(m)]);
    
    xdat        = modelz{m};
    
    for g = 1:ngroups
        
        % Omit first trial (doesn't have subsequent update, so can't be estimated)
        indx1 = ones(ttrials,1);
        indx1(1:btrials:ttrials) = 0;
        
        idx = data.group == groups(g) & indx;
        
        % Regression
        [out, stats, mdl] = regr_fitglm(xdat(idx,:),ydat(idx),data.sub(idx),rtype,ptype,0,normalise,[],indx1(idx));
        
        % Save LL per group
        if g == 1
            LL_patient(:,m) = out.LL;
        else
            LL_control(:,m) = out.LL;
        end
        
        % Get number of regressors used
        nregr   = size(out.Betas,2)-1;
    end
end

bms_conf(:,1) = spm_BMS(LL_patient);
bms_conf(:,2) = spm_BMS(LL_control);

% Plot BMS results
if do.plotting
    
    % Plot
    figC = figure; hold on;
    bhand = bar(bms_conf);
    bhand(1).FaceColor = 'r';
    bhand(1).EdgeColor = 'none';
    bhand(2).FaceColor = [1 1 1]*.9;
    legend(bhand,grouplab,'FontSize',titlefntsz,'Location','NorthWest');
    legend boxoff
    set(gca,'XTick',1:4,'XTickLabel',{'Model 1','Model 2','Model 3','Model 4'},'FontSize',axfntsz);
    yLab = ylabel('Posterior model probabilities'); set(yLab,'FontSize',titlefntsz);
    
    if do.saving_plot
        save2eps(figC,figpath,'BMS_confidence',[.1 .1 600 400]);
    end
        
end

%% Run different regression models: update

modelz  = {[X1],[X1 X1.*X2],[X1 X1.*X2 X1.*(1-X2).*X3],[X1 X1.*X2 X1.*(1-X2).*X3 X1.*X4]};

ydat        = data.lr.*abs(data.delta);
rtype       = 'normal';
ptype       = 'linear';
normalise   = 1;

regrlab = {'|PE|','CPP','(1-CPP)*(1-MC)','hit'};

for m = 1:size(modelz,2)
    
    disp(['Running regression model ' num2str(m)]);
    
    xdat        = modelz{m};
    
    for g = 1:ngroups
        
        % Omit first trial (doesn't have subsequent update, so can't be estimated)
        indx1 = ones(ttrials,1);
        indx1(1:btrials:ttrials) = 0;
        
        idx = data.group == groups(g) & indx;
        
        % Regression
        [out, stats, mdl] = regr_fitglm(xdat(idx,:),ydat(idx),data.sub(idx),rtype,ptype,0,normalise,[],indx1(idx));
        
        % Save LL per group
        if g == 1
            LL_patient(:,m) = out.LL;
        else
            LL_control(:,m) = out.LL;
        end
        
        % Get number of regressors used
        nregr   = size(out.Betas,2)-1;
    end
end

% Calculate BMS (iterative)
bms_update(:,1) = spm_BMS(LL_patient);
bms_update(:,2) = spm_BMS(LL_control);

% Plot different BIC
if do.plotting
    
    % Plot
    figU = figure; hold on;
    bhand = bar(bms_update);
    bhand(1).FaceColor = 'r';
    bhand(1).EdgeColor = 'none';
    bhand(2).FaceColor = [1 1 1]*.9;
    legend(bhand,grouplab,'FontSize',titlefntsz,'Location','NorthWest');
    legend boxoff
    set(gca,'XTick',1:4,'XTickLabel',{'Model 1','Model 2','Model 3','Model 4'},'FontSize',axfntsz);
    yLab = ylabel('Posterior model probabilities'); set(yLab,'FontSize',titlefntsz);
      
    if do.saving_plot
        save2eps(figU,figpath,'BMS_update',[.1 .1 600 400]);
    end
    
end