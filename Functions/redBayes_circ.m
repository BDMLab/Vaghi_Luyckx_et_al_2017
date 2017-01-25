function [B,lr,d,CPP,MC,r,sF] = redBayes_circ(x,H,sN,low,up)
% function [B,lr,d,CPP,MC,r,sF] = redBayes_circ(x,H,sN,low,up)
%   
% Input:
%   - x     = vector with values to estimate mean
%   - H     = hazard rate
%   - sN    = std of normal distribution (generative std)
%   - low   = lower bound of uniform distribution
%   - up    = upper bound of uniform distribution
%
% Output:
%   - B     = beliefs about underlying distribution (= estimate)
%   - lr    = learning rate
%   - d     = prediction error between current belief and observation
%   - CPP   = change point probability
%   - MC    = model confidence (= 1-RU)
%   - r     = expected run length
%   - sF    = expected amount of noise from the gen distribution
%
% Fabrice Luyckx, 20/3/2016

% Get dimensions
sz      = size(x);
ttrials = sz(1);
nsamp   = sz(2);

% Initialise variables
B       = zeros(ttrials,nsamp); % belief about current mean
lr      = B; % learning rate
d       = B; % prediction error
CPP     = B; % change-point probability (omega)
MC      = B+0.5; % relative uncertainty (tau)
r       = B+1; % expected run length
sF      = B; % expected amount of noise from the gen distribution

for t = 1:ttrials
    for i = 1:nsamp
        
        % Prediction error (observed point - belief)
        d(t,i)      = diffcirc(B(t,i),x(t,i),1);
                
        % Change-point probability
        %sF(t,i)     = sqrt(sN^2 + (sN^2)/r(t,i)); %2010
        sF(t,i)     = sqrt(sN^2 + (sN^2 * (1-MC(t,i)))/(MC(t,i))); %2012
        CPP(t,i)	= (unifpdf(x(t,i),low,up)*H)/((unifpdf(x(t,i),low,up)*H) + (normpdf(diffcirc(x(t,i),B(t,i)),0,sF(t,i))*(1-H)));
        
        % Learning rate
        lr(t,i)     = CPP(t,i) + (1-MC(t,i))*(1-CPP(t,i));
        
        % Next trial calculations
        
        if i < nsamp
            % Belief on next trial
            B(t,i+1)  = mod(B(t,i) + lr(t,i)*d(t,i),360);
            
            % MC on next trial
            term1   = CPP(t,i)*sN^2; % weighted average of the variance on the distribution conditional on a change point
            term2   = (1-CPP(t,i))*(1-MC(t,i))*sN^2; % weighted average of the variance on the distribution conditional on NO change point
            term3   = CPP(t,i)*(1-CPP(t,i))*(d(t,i)*(MC(t,i)))^2; % accounts for variance emerging from the difference in the means of these two conditional distributions
            MC(t,i+1) = 1 - ((term1+term2+term3)/(term1+term2+term3 + sN^2)); % denominator contains an additional term to account for uncertainty arising from noise
            
            % Expected run length on next trial
            r(t,i+1)  = (r(t,i) + 1)*(1-CPP(t,i)) + CPP(t,i);
        end        
    end
end

end

%% Differencing degrees

function [difference] = diffcirc(a,b,varargin)
%function [difference] = diffcirc(a,b)
%
% Function to calculate difference between two angles (in degrees).
% Optional: 'direction', if 1 gives the sign of the difference (a > b?)
%
% Fabrice Luyckx, 8/4/2016

% Optional arguments
optargs = {0};
specif = find(~cellfun(@isempty,varargin));
[optargs{specif}] = varargin{specif};
[direction] = optargs{:};

% Find smallest and largest value
c       = min([a,b]);
d       = max([a,b]);

% Possible distances
opt1    = mod(c - d,360);
opt2    = mod(d - c,360);

% Shortest difference between angles
difference = min([opt1 opt2]);

% Difference (with or without sign)
if direction == 1
    if abs(a-b) < 180
        difference = -1*sign(a-b)*difference;
    else
        difference = sign(a-b)*difference;
    end
end

end
