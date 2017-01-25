function [lr, delta] = estimateLR(samples,bposition)
% function [lr, delta] = estimateLR(samples,bposition)
%
% Estimate learning rate and delta from data.
%
% In:
%   - samples: (block x blocktrials)
%   - bposition: (block x blocktrials)
%
% Out: 
%   - lr: estimated learning rate
%   - delta: estimated prediction error
%
% Fabrice Luyckx, 27/9/2016

%% Get some values

lr      = 0*bposition;
delta  	= 0*bposition;

rowz    = size(bposition,1);
columnz = size(bposition,2);

%% Estimate learning rates

for b = 1:rowz
    for t = 1:columnz
        if t < columnz
            lr(b,t)    = (diffcirc(bposition(b,t),bposition(b,t+1)))/(diffcirc(bposition(b,t),samples(b,t)));
            % When coin lands exactly in middle of bucket
            if abs(lr(b,t)) == Inf
                lr(b,t)    = diffcirc(bposition(b,t),bposition(b,t+1));
            end
        end
        delta(b,t) = diffcirc(bposition(b,t),samples(b,t));
    end
end

lr(isnan(lr)) = 0;

end

