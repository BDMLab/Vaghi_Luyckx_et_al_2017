%%%%%%%%%%%%%%%%
%% Instructions
%%%%%%%%%%%%%%%%

switch instructions
    
    %% Main instructions
    case 'beforepractice'
        
        % Welcome screen
        Screen(w,'FillRect',col.background);   % background colour
        Screen('TextStyle',w, 1); % bold
        Screen('TextSize', w , titleSize);
        Screen(w,'DrawText','Welcome!',leftMargin,screenYpixels * 0.1,col.white);
        
        Screen('TextSize', w,textSize);
        Screen('TextStyle',w, 0); %normal
        Screen(w,'DrawText','Press spacebar to continue ...',leftMargin,screenYpixels * 0.9,col.white);
        
        Screen(w,'Flip');
        
        KbWait;
        while KbCheck; end;
        
        % //Page 1
        
        Screen(w,'FillRect',col.background);
        mytext = printText('Introduction_1.m',10);        
        DrawFormattedText(w, mytext, leftMargin, topMargin, col.white,[],[],[],1.5);       
        Screen(w,'DrawText','Press spacebar to continue ...',leftMargin,screenYpixels * 0.9,col.white);       
        Screen(w,'Flip');
        
        KbWait;
        while KbCheck; end;
        
        % //Page 2
        mytext = printText('Introduction_2.m',13);        
        DrawFormattedText(w, mytext, leftMargin, topMargin, col.white,[],[],[],1.5);

        %% Instructions after practice trials
    case 'afterpractice'
                
        Screen(w,'FillRect',col.background);
        Screen('TextSize', w,textSize);
        Screen('TextStyle',w,0); % normal
        
        mytext = printText('Conclusion_1.m',13);
        DrawFormattedText(w, mytext, leftMargin, topMargin, col.white,[],[],[],1.5);
        
        %% Instructions at start of each block
    case 'startblock'
        
        Screen(w,'FillRect',col.background);
        Screen('TextSize', w,titleSize);
        Screen('TextStyle',w,1); % bold
        
        if practice == 0
            line1 = ['Run ' num2str(data.block(t))];
            DrawFormattedText(w, line1,leftMargin, screenYpixels * 0.1, col.white,[],[],[],1.5);
        else
            line1 = 'Practice run';
            DrawFormattedText(w, line1,leftMargin, screenYpixels * 0.1, col.white,[],[],[],1.5);            
        end
        
        %% Instructions at the end of each block
    case 'endblock'
        
        Screen(w,'FillRect',col.background);
        Screen('TextSize', w,titleSize);
        Screen('TextStyle',w,1); % bold
        line1 = 'Run completed';
        DrawFormattedText(w, line1,leftMargin, screenYpixels * 0.1, col.white,[],[],[],1.5);
        
        Screen('TextSize', w,textSize);
        Screen('TextStyle',w, 0); % normal
        line1 = 'Take a short break.';
        line2 = ['\nTotal points earned so far: ' num2str(round(data.totReward(t)))];
        DrawFormattedText(w,[line1 line2], 'center', 'center', col.white,[],[],[],1.5);
        
        %% Instructions at end of experiment
    case 'endexp'
        
        Screen(w,'FillRect',col.background);
        Screen('TextSize', w,titleSize);
        Screen('TextStyle',w,1); % bold
        line1 = 'Experiment finished!';
        DrawFormattedText(w, line1,'center', screenYpixels * 0.25, col.white,[],[],[],1.5);
        
        Screen('TextSize', w,textSize);
        Screen('TextStyle',w, 0); % normal
        line1 = ['You earned a total of ' num2str(round(data.totReward(t))) ' points'];
        line2 = ['\nYou''ve earned ' num2str(money) ' pounds.'];
        DrawFormattedText(w,[line1 line2], 'center', 'center', col.white,[],[],[],1.5);
        
end

Screen('TextSize', w,textSize);
Screen('TextStyle',w, 0); % normal
Screen(w,'DrawText','Press spacebar to continue ...',leftMargin,screenYpixels * 0.9,col.white);
Screen(w,'Flip');

KbWait;
while KbCheck; end;

Screen(w,'FillRect',col.background);
Screen(w,'Flip');
WaitSecs(.5);