function [allAngles] = circle(numslots, circlerad,plotornot)
% function [allAngles] = circle(numslots, circlerad,plotornot)
%
% Function to get n angles from circle with certain radius.
%
% Input:
%   numslots    = number of angles
%   circlerad   = radius [200]
%   plotornot   = plot it? [0]
%
% Output:
%   allAngles = angles

if nargin<2
    circlerad = 200;
elseif nargin<3
    plotornot = 0;
end

theta = 360/numslots;

if plotornot == 1
    figure;
end

for n = 1:numslots;
    
    beta    = 90-(n*theta);
    x       = cos((beta+theta)*pi/180).*circlerad;
    y       = sin((beta+theta)*pi/180).*circlerad;
    a       = cos((beta)*pi/180).*circlerad;
    b       = sin((beta)*pi/180).*circlerad;
    
    for c = 1:length(circlerad)
        angles(c,n,:)=[x(1) y(1) a(1) b(1)];
    end
      
    if plotornot == 1
        hold on;
        for c = 1:length(circlerad)
            line([angles(c,n,1) angles(c,n,3)],[angles(c,n,2) angles(c,n,4)])
        end
        
        plot(0,0,'o');
        set(gca,'Xlim',[-max(circlerad)*2 max(circlerad)*2])
        set(gca,'Ylim',[-max(circlerad)*2 max(circlerad)*2])
        axis square
    end
    
    allAngles = squeeze(angles(:,:,1:2));
    
end