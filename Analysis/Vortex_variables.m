%%%%%%%%%%%%%%%%%%%
%% DEFINE VARIABLES
%%%%%%%%%%%%%%%%%%%

% Generates general variables needed in subsequent scripts. Defines plot
% variables that might be used in later plots.

submat      = unique(data.sub)';
nsubj       = length(submat);
gsubj       = [length(unique(data.sub(data.group == 0))), length(unique(data.sub(data.group == 1)))]; % nsubj per group
ttrials     = length(data.sub);
ntrials     = length(data.sub(data.sub == submat(1)));
nblocks     = length(unique(data.block));
btrials     = ntrials/nblocks;
groups      = unique(data.group)';
ngroups     = length(groups);

part.group  = str2num(participant.group);

gsubmat{1} = unique(data.sub(data.group == 0))';
gsubmat{2} = unique(data.sub(data.group == 1))';

%% Functions

% Function to make matrix one long column
makeLong = @(x) x(:);

% Function to calculate color matrices
rgb     = @(x) round(x./255,2);

%% Plot variables

axfntsz         = 14;
titlefntsz      = 16;
lnwid           = 2;
mksz            = 8;
barwid          = .9;

set(0,'DefaultAxesFontName', 'Helvetica');
set(0,'DefaultTextFontname', 'Helvetica');

modmark     = {'d','x'};
markers 	= {'o','s','d','x','h','p','+','*'};

rbocol  = rgb([101, 153, 255]);
humcol  = rgb([255, 153, 0]);
CPPcol 	= rgb([101, 153, 255]);
RUcol  	= [.2 .4 .53];
lrcol   = rgb([175, 234, 170]);

col = jet(nsubj);

% colz(1,:)   = [0.2, 0.8, 0.5];
% colz(2,:) 	= [0, 0.5, .75];
% colzedge    = colz.*.8;

colz(1,:)   = rgb([183, 110, 184]);
colz(2,:) 	= rgb([122, 186, 122]);
colzedge    = colz.*.8;


grouplab    = {'Patients','Controls'};
