%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vortex: hazard rate fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%% Run model with different hazard rates

% Belief updates (human) / belief updates model = lr*d
samples             = reshape(stim.samples,btrials,nblocks*nsubj)'; % observed data
bposition           = reshape(data.bucketPosition,btrials,nblocks*nsubj)'; % bucket position
submat2             = makeLong(repmat(submat',1,nblocks)');

% Model input
x      	= reshape(stim.samples,btrials,nblocks*nsubj)'; % observed data
sN      = stim.std(1); % standard deviation normal distribution
low     = 1;
up      = stim.nPositions(1);

Hspace  = [0:.01:1]; % hazard rate
B       = zeros(nblocks*nsubj,btrials,length(Hspace));
CPP     = B;

% Reduced Bayesian obsession
for m = 1:length(Hspace)
    
    disp(['H = ' num2str(Hspace(m))]);
    
    [B(:,:,m),~,~,CPP(:,:,m)] = redBayes_circ(x,Hspace(m),sN,low,up);
    
    for s = 1:nsubj
        
        idx = submat2 == submat(s);
        
        Bmod(:,m,s) = makeLong(B(idx,2:end,m)');
        bhum(:,s)	= makeLong(bposition(idx,2:end)');
        
        for t = 1:length(Bmod)
            sqDiff(t) = diffcirc(Bmod(t,m,s),bhum(t,s)).^2;
        end
        
        Lsq(m,s) = sum(sqDiff);
    end
end

[~,minLsq]  = min(Lsq);
Hfits        = Hspace(minLsq);

%% Test significance

if do.signif    
    gidx = str2num(participant.group);  
    [p,h] = ranksum(Hfits(gidx==1),Hfits(gidx==0));
    median(Hfits(gidx==1))
    median(Hfits(gidx==0))
end

%% PLOTS

if do.plotting
    
    %% Distribution of all fits
    figure; histogram(Hspace(minLsq),length(Hspace));
    
    %% Individual fits over whole space
    figure;
    plot(Hspace,Lsq); xlim([Hspace(1) Hspace(end)]);
    ylabel('Sum of squares');
    xlabel('Parameter space');
    title('LS for each subject');
    
    %% Distribution of two groups
    
    figD = figure;
    nbins = 10;
    gidx = str2num(participant.group);
    
    subplot(1,2,1); hold on;
    h(1) = histogram(Hfits(gidx==1),nbins,'FaceColor',[1 1 1].*.9,'FaceAlpha',1);
    title('Controls');
    ylim([0 10]);
    xlim([0 1])
    ylabel('n subjects','FontSize',axfntsz);
    xlabel('Parameter space','FontSize',axfntsz);
    
    subplot(1,2,2); hold on;
    h(2) = histogram(Hfits(gidx==0),nbins,'FaceColor','r','FaceAlpha',1);
    title('Patients');
    ylim([0 10]);
    xlim([0 1])
    ylabel('n subjects','FontSize',axfntsz);
    xlabel('Parameter space','FontSize',axfntsz);
    %title('Distribution of parameter fits','FontSize',14);
    
    if do.saving_plot
        save2eps(figD,figpath,'Hazard_fits',[.1 .1 600 400])
    end
    
end

%% Estimate learning rate 

% //REDUCED BAYESIAN OBSERVER

B           = zeros(nsubj*nblocks,btrials);
lr          = B;
d           = B;
MC          = B;
CPP         = B;
r           = B;

for s = 1:nsubj
    
    idx = submat2 == submat(s);
    
    H = Hfits(s); % hazard rate
    [B(idx,:),lr(idx,:),d(idx,:),CPP(idx,:),MC(idx,:),r(idx,:)] = redBayes_circ(x(idx,:),H,sN,low,up);
end

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
X5      = data.zDiffBPos'; %absolute difference in reported confidence 

shiftX1 = shift(X1,1,data.sub,data.block);
shiftX2 = shift(X2,1,data.sub,data.block);
shiftX4 = shift(X4,1,data.sub,data.block);

%% Regression: confidence

xdat        = [shiftX1 shiftX2 (1-shiftX2).*X3 shiftX4];

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
    
    %Save stats
    reg_RegCF{g}     = stats;
    out_RegCF{g}     = out;
    
    %Build file for Betas to be plotted
    Beta    = [Beta;out.Betas];
    G       = [G;repmat(groups(g),size(out.Betas,1),1)];
    
end

% Plot

if do.plotting
        
    ptnt            = Beta(G==0,2:end);
    exclude_ptnt    = [1,4,15,19];
    valid_ptnt      = setdiff(1:24,exclude_ptnt);
    m_ptnt          = ptnt(valid_ptnt,:);
    
    sem_ctrl        = std(Beta(G==1,2:end))./sqrt(gsubj(1));
    sem_ptnt        = std(m_ptnt)./sqrt(gsubj(2)-4);
    
    figC = figure; hold on;
    plot([0 5],[0 0],'k--');
    cplot = errorbar(1:4,mean(Beta(G==1,2:end)),sem_ctrl,'o',...
        'LineWidth',lnwid,'MarkerFaceColor',[1 1 1].*.8,'MarkerEdgeColor',[1 1 1].*.6,'Color',[1 1 1].*.6,'MarkerSize',10);
    pplot = errorbar(1:4,mean(ptnt(valid_ptnt,:)),sem_ptnt,'o',...
        'LineWidth',lnwid,'MarkerFaceColor','r','MarkerEdgeColor',[1 0 0].*.85,'Color',[1 0 0].*.85,'MarkerSize',10);
    legend([cplot,pplot],{'Controls','Patients'},'Location','NorthWest','FontSize',titlefntsz);
    legend boxoff
    set(gca,'XTick',1:4,'XTickLabel',regrlab,'FontSize',axfntsz);
    
    xlabel('Predictors','FontSize',axfntsz);
    ylabel({'Regression coefficients';'beta weights'},'FontSize',axfntsz);
    
    %title('Confidence regression');
    
    if do.signif
        [h,p,ci,stats] = ttest2(Beta(G==1,2:end),ptnt(valid_ptnt,:));
        h(h==0) = nan;
        
        plot(1:4,h*.1,'k*','MarkerSize',10,'LineWidth',2);        
    end
    
    if do.saving_plot
        save2eps(figC,figpath,'Confidence_regressions_Hfit',[.1 .1 600 400])
    end
    
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
    
    %Save stats
    reg_RegLR{g}    = stats;
    out_RegLR{g} 	= out;
    
    %Build file for Betas to be plotted
    Beta            =[Beta;out.Betas];
    G               =[G;repmat(groups(g),size(out.Betas,1),1)];
end

% Plot

if do.plotting
        
    ptnt            = Beta(G==0,2:end);
    exclude_ptnt    = [1,4,15,19];
    valid_ptnt      = setdiff(1:24,exclude_ptnt);
    m_ptnt          = ptnt(valid_ptnt,:);
    
    sem_ctrl        = std(Beta(G==1,2:end))./sqrt(gsubj(1));
    sem_ptnt        = std(m_ptnt)./sqrt(gsubj(2)-4);
    
    figL = figure; hold on;
    plot([0 5],[0 0],'k--');
    cplot = errorbar(1:4,mean(Beta(G==1,2:end)),sem_ctrl,'ko','MarkerFaceColor',[1 1 1].*.9);
    pplot = errorbar(1:4,mean(ptnt(valid_ptnt,:)),sem_ptnt,'ko','MarkerFaceColor','r');
    legend([cplot,pplot],{'Controls','Patients'});
    set(gca,'XTick',1:4,'XTickLabel',regrlab,'FontSize',axfntsz);
    title('Update regression');
end