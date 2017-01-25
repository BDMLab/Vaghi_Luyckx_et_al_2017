function [out, stats, mdl] = regr_fitglm(xdat,ydat,subz,varargin)
% function [out, stats, mdl] = regr_fitglm(xdat,ydat,subz,[rtype, ptype, plotting, normalise, condz, indx, intercept])
%
% Function to run regressions per subject.
%   
% Input:
%   xdat        = regressors
%   ydat        = dependent variable
%   subz        = subject indices   (= vector length of whole experiment)
%   [rtype]     = regression type ['normal']
%   [ptype]     = values type: e.g. constant, linear, interactions, ... ['linear']
%   [plotting]  = plot or not? [0];
%   [normalise] = normalise predictors? [1]
%   [condz]     = vector indication conditions 1 - n [ones(ttrials,1)]
%   [indx]      = subset, for example to exclude bad trials [ones(ttrials,1)]
%   [intercept] = plot with intercept? [0]
%
% Output:
%   out:
%       - Betas     = Beta coefficients
%       - Err       = standard error of betas
%       - LL        = loglikelihood of each regression
%       - rsquared  = r square of each regression
%
%   stats:
%       - h         = reject null-hypothesis
%       - p         = p-values for each Beta weight
%
%   mdl: example model output from last subject (and condition)
%
%  Fabrice Luyckx, 13/4/2016

%% Get some values

submat  = unique(subz)';
nsubj   = length(submat);
ttrials = length(subz);

%% DEFAULT VALUES

optargs = {'normal' 'linear' 0 1 ones(ttrials,1) ones(ttrials,1) 0};

% Now put these defaults into the valuesToUse cell array,
% and overwrite the ones specified in varargin.
specif = find(~cellfun(@isempty,varargin)); % find position of specified arguments
[optargs{specif}] = varargin{specif};

% Place optional args in memorable variable names
[rtype, ptype, plotting, normalise, condz, indx, intercept] = optargs{:};

ncond = length(unique(condz));

fail = []; % variable to store failed fits
warning('');

%% NORMALISE the regressors

% Normalise per subject
if normalise == 1
    for s = 1:nsubj
        xdat(subz == submat(s),:) = zscore(xdat(subz == submat(s),:));
    end
end

%% REGRESSION

disp(' ');
disp('Running regression ...');

for c = 1:ncond
    
    if ncond > 1
        disp(' ');
        disp(['Condition ' num2str(c)]);
    end
    disp(' ');
    
    % REGRESSION
    for s = 1:nsubj
        disp(num2str(s));
        index = find(subz == submat(s) & condz == c & indx);
               
        mdl                 = fitglm(xdat(index,:),ydat(index),ptype,'Distribution',rtype);
        
        out.Betas(s,:,c)    = mdl.Coefficients.Estimate;
        out.LL(s,c)         = mdl.LogLikelihood;
        out.rsquared(s,c)   = mdl.Rsquared.Adjusted;
        %out.residuals(s,:,c) = mdl.Residuals.Raw;
        
        % Catch warning
        msgstr = lastwarn;
        if ~isempty(msgstr)
            fail = [fail, s];
            warning('');
        end
    end
    
    exclude                 = sort(unique(fail));
    out.Betas(exclude,:,:)  = nan;
    if ~isempty(exclude)
        warning(['Excluded participants: ' num2str(exclude)]);
    end
    
    % standard error
    out.Err = nanstd(out.Betas)/sqrt(nsubj);   
end

%% STATS
[stats.h, stats.p, stats.ci, stats.stats] = ttest(out.Betas);

%% PLOT

if plotting == 1
    
    fig = gcf; hold on;

    % Colours
    if ncond == 1
        colz = 'k';
    else
        colmap  = winter;
        colz    = colmap(linspace(1,length(colmap),ncond),:);
    end
        
    % Plot with intercept?
    if intercept == 1
        pB  = out.Betas;
        pE  = out.Err;
    else
        pB  = out.Betas(:,2:end,:);
        pE  = out.Err(:,2:end,:);
    end
    
    % Actual plot
    for c = 1:ncond
        hand(c) = errorbar(nanmean(pB(:,:,c),1),1.96.*pE(1,:,c),'o','Color',colz(c,:),'LineWidth',2,'MarkerSize',7,'MarkerFaceColor',colz(c,:)); % neutral
    end
   
    % Plot 0 line
    plot([0 size(xdat,2)+1],[0 0],'k--');
    
end

end

