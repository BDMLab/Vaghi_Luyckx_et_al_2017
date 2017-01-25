function [out] = shift(xdat,n, subz, blockz, varargin)
%function [out] = shift(xdat,n)
%
% Function to shift vector with n places further.
%
% - Input:
%   xdat        = matrix to shift (columns)
%   n           = number of places to shift
%   subz        = subject index
%   block       = block index
%   [direction] = 'columns' or 'rows'
%
% - Output:
%   out     = shifted vector
%
% Fabrice Luyckx, 1/04/2016

if isempty(varargin)
    sz = size(xdat);
    if sz(2) == 1
        direction = 'columns';
    elseif sz(1) == 1
        direction = 'rows';
    end
else
    direction = varargin{1};
end

submat  = unique(subz);
nsubj   = length(submat);
bmat    = unique(blockz);
nblocks = length(bmat);

% Make shift
out             = zeros(size(xdat));

for s = 1:nsubj
    for b = 1:nblocks
        idx = subz == submat(s) & blockz == bmat(b);
        tmp = xdat(idx);
        
        switch direction
            case 'columns'
                out(idx,:) = [zeros(n,1); tmp(1:end-n,:)];
            case 'rows'
                out(:,idx) = [zeros(1,n) tmp(:,1:end-n)];
        end
    end
end

end

