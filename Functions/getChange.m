function [changePatterns, cErr, changePoint, subz2] = getChange(dat,changez,lag,subz,varargin)
% function [changePatterns, cErr, changePoint] = getChange(dat,changez,lag,subz,[useMean],[condz],[indx])
%
% Function to find change points in data
%
% Input:
%   - dat           = data to average over
%   - changez       = zero vector with changes (1)
%   - lag           = how many positions before and after the change point
%   - subz          = subject indices
%   - [useMean]     = use mean (1) or median (0) summary [1]
%   - [condz]       = conditions [ones(ntrials,1)]
%   - [indx]        = subset data [ones(ntrials,1)]
%
% Output:
%   - changePatterns    = mean of dat on each position
%   - cErr              = standard error for each position
%   - changePoint       = change point matrix
%   - subz2             = subject index for change point matrix
%
% Fabrice Luyckx, 9/5/2016

%% Get some values

submat  = unique(subz)';
nsubj   = length(submat);
ttrials = length(subz);

%% DEFAULT VALUES

optargs = {1,ones(ttrials,1) ones(ttrials,1)};

% Now put these defaults into the valuesToUse cell array,
% and overwrite the ones specified in varargin.
specif = find(~cellfun(@isempty,varargin)); % find position of specified arguments
[optargs{specif}] = varargin{specif};

% Place optional args in memorable variable names
[useMean, condz, indx] = optargs{:};

cmat    = unique(condz);
ncond   = length(cmat);

%%

% Find points where switch has occured
pattern         = [zeros(1,lag) 1 zeros(1,lag)]; % pattern to look for
patternIndices  = strfind(changez,pattern);
excludePattern  = [];

% Exclude patterns with outlying value in it (indx)
if sum(find(indx==0)) >= 1
    
    badTrials       = find(indx==0);
    fullPatternInd  = zeros(length(patternIndices),length(pattern));
    
    for i = 1:length(patternIndices)
        fullPatternInd(i,:) = patternIndices(i):patternIndices(i)+length(pattern)-1;
        
        for b = 1:length(badTrials)
            if any(fullPatternInd(i,:) == badTrials(b))
                excludePattern = [excludePattern, i];
                break;
            end
        end
    end   
end

tmp = patternIndices(setdiff(1:length(patternIndices),excludePattern));

changePoint     = zeros(length(tmp),length(pattern));
subz2           = zeros(length(tmp),1);
btypez          = zeros(length(tmp),1);
changePatterns  = zeros(nsubj,length(pattern),ncond);

% Track data for change points
for i = 1:length(tmp)
    changePoint(i,:) 	= dat(tmp(i):tmp(i)+length(pattern)-1);
    subz2(i)            = subz(tmp(i));
    btypez(i)         	= condz(tmp(i));
end

% Get pattern per subject and condition
for s = 1:nsubj
    for c = 1:ncond
        idx = subz2 == submat(s) & btypez == cmat(c);  
        if useMean == 0
            changePatterns(s,:,c) = nanmedian(changePoint(idx,:));
        else
            changePatterns(s,:,c) = nanmean(changePoint(idx,:));
        end
    end
end

% SEM
cErr = mad(changePatterns);

end

