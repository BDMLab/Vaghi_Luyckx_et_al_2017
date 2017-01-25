%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vortex: error magnitude
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generates data presented in Figure 1D and writes it to a data file that
% can produce the plot in R.

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
do.writeR           = true; % write data to R file for plots?

%% Set variables

% //REDUCED BAYESIAN OBSERVER

x      	= reshape(stim.samples,btrials,nblocks*nsubj)'; % observed data
H       = stim.changeProb(1); % hazard rate
sN      = stim.std(1); % standard deviation normal distribution
low     = 1;
up      = stim.nPositions(1);

% Belief updates (human) / belief updates model = lr*d
samples             = reshape(stim.samples,btrials,nblocks*nsubj)'; % observed data
bposition           = reshape(data.bucketPosition,btrials,nblocks*nsubj)'; % bucket position

[data.lr, data.delta] = estimateLR(samples,bposition);

% Reshape for plotting
data.lr         = makeLong(data.lr');
data.delta      = makeLong(data.delta');
data.confidence = makeLong(data.confidence');

% //SET INDEX
% Exclude trials where estimated learning rate is higher than 99% of group
% estimates (extreme outliers)and where delta is equal 0 (see Nassar 2016)
indx    = [data.lr(data.group==1)<prctile(data.lr(data.group==1),99); data.lr(data.group == 0)<prctile(data.lr(data.group==0),99)]& data.delta > 0;

%% Bin errors/learning rate

nbins       = 20; %to have quant space based on data.delta excluding 0
tmp         = data.delta(data.delta~=0);
quantSpace  = [min(tmp) quantile(tmp,nbins-1) max(tmp)];

xValues     = quantSpace(1:end-1)+diff(quantSpace)/2;

for s = 1:nsubj
    for q = 1:nbins
        idx         = data.sub == submat(s) & data.delta >= quantSpace(q) & data.delta < quantSpace(q+1) & indx;
        qLr(s,q)    = mean(data.lr(idx));
    end
end

for g = 1:ngroups
    lrMeans(g,:)    = mean(qLr(part.group == g-1,:));
    %this for CI
    % lrSem(g,:)      = 1.96*std(qLr(part.group == g-1,:))./sqrt(sum(qLr(part.group == g-1,:) ~= 0));
    lrSem(g,:)      = std(qLr(part.group == g-1,:))./sqrt(sum(qLr(part.group == g-1,:) ~= 0));
end

%Bins
ErrBins=[(1:1:20)'; (1:1:20)']
%Means per group based on 20 bins
ErrMeans= [lrMeans(1,:)'; lrMeans(2,:)'];
%SEM per group based on 20 bins
ErrSem=   [lrSem(1,:)'; lrSem(2,:)'];
% Code groups
ErrGroup=[repmat(0,20,1); repmat(1,20,1)];
%X axis
ErrxValues= [xValues'; xValues'];

%% Reshape data for plots in R

if do.writeR
    nbins       = 3;
    %to have quant space based on data.delta excluding 0
    tmp         = data.delta(data.delta~=0);
    quantSpace  = [min(tmp) quantile(tmp,nbins-1) max(tmp)];
    
    [ind1,~]=find(ErrxValues(:,:)>quantSpace(1) & ErrxValues(:,:)<=quantSpace(2));
    [ind2,~]=find(ErrxValues(:,:)>quantSpace(2) & ErrxValues(:,:)<=quantSpace(3));
    [ind3,~]=find(ErrxValues(:,:)>quantSpace(3) & ErrxValues(:,:)<=quantSpace(4));
    
    ErrMagn=[ErrBins, ErrxValues, ErrMeans, ErrSem, ErrGroup];
    ErrMagn(ind1,6)=1;
    ErrMagn(ind2,6)=2;
    ErrMagn(ind3,6)=3;
    
    filename=fullfile(datapath, 'RPlotData', 'VortexBet_ErrorMag.txt')
    dlmwrite(filename, ErrMagn)
end


%% based on 3 bins, low, high , medium 

keySet      = submat;
valueSet    = str2num(participant.group);
mapObj      = containers.Map(keySet, valueSet);


nbins       = 3;
%to have quant space based on data.delta excluding 0
tmp         = data.delta(data.delta~=0);
quantSpace  = [min(tmp) quantile(tmp,nbins-1) max(tmp)];

for s = 1:nsubj    
    for q = 1:nbins
        idx         = data.sub == submat(s) & data.delta >= quantSpace(q) & data.delta < quantSpace(q+1) & indx;
        qLr(s,q)    = mean(data.lr(idx));
            qLr(s,4)    = submat(s)
       qLr(s,5)    =  mapObj(submat(s));
 
    end   
end

if do.writeR
    filename = fullfile(datapath, 'RPlotData', 'ErrorMagn_3bins.txt');
    dlmwrite(filename, qLr);
end
   