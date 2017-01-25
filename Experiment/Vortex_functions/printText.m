function [mytext] = printText(textfile,finishLine)
% function [mytext] = printText(textfile,finishLine)
%   
% Function to read out and print text from a textfile with PTB.

fd = fopen(textfile, 'rt');
if fd==-1
    error('Could not open file!');
end

mytext = '';
tl = fgets(fd);
lcount = 0;
while lcount < finishLine
    mytext = [mytext tl];
    tl = fgets(fd);
    lcount = lcount + 1;
end
fclose(fd);

% Get rid of '% ' symbols at the start of each line:
mytext = strrep(mytext, '% ', '');
mytext = strrep(mytext, '%', '');

end

