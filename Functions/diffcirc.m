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

