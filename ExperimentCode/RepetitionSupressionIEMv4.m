% function RepetitionSupressionIEMv4

try
    clc; clear all; close all;
    addpath('C:\Users\PsychEEGlab\Desktop\EXPERIMENTS\Matt\Repetition supression')
    
    %% set up the save file names
    SubNo= input('Enter the subject number \n','s');

    CountBalance=input('Counter-balance order (1-2) \n','s');
    CountBalance=str2double(CountBalance);
    now=clock; % gets the current time and date for unique file name
    dataFileDefault=['Sub',SubNo,'(', num2str(now(3)),'-', num2str(now(2)),'-',...
        num2str(now(4)),num2str(now(5)),')','-RepetitionSupressionIEMv4'];
    dataFileName=strcat(dataFileDefault,'.csv');
    dataFileNameBU=strcat(dataFileDefault,'-BU.csv');
    
    formatString = repmat('%d, ',1,18);
    formatString=strcat(formatString,' \n');
    SubNo=str2double(SubNo);

    RedMod=.7;
    GreenMod=.7;
    input(' \n ********************************\n    Is the EEG recording data?\n ********************************')
    
    %% set up PT
    PsychDefaultSetup(2);
    bg=0.5; %mid-gray background luminance
    Screen('Preference', 'SkipSyncTests', 0);
    screenNum=2;
    [w ,rect]=PsychImaging('OpenWindow',screenNum,bg);
    ifi=Screen('GetFlipInterval', w);
    Hz=120; % refresh rate in Hz
    center=[rect(3)/2  rect(4)/2];
    Priority(MaxPriority(w));
    %     HideCursor;
    KbName('UnifyKeyNames'); % makes the keyboard compatible between systems
    rng shuffle; % shuffles the random number generator
    %% sets up the EEG triggers
    pulse_duration = 4; % This will be how long the trigger is sent for, in milliseconds.
    
    parallel_port_address = hex2dec('D050');
    send_triggers =1;
    % send a trigger to a non-existent parallell port.
    if send_triggers
        ioObj = io64; % create an instance of the io64 object
        status = io64(ioObj); % initialize the inpout64 system driver
    end
    if send_triggers
        io64(ioObj,parallel_port_address,0); % set the port called parallel_port_address to 0
    end
    %% makes sure the monitor is correctly set up
    if send_triggers
        if round(1/ifi)~=Hz
            sca
            clc
            disp('wrong refresh rate')
        end
        if rect(3)~=1920 || rect(4) ~=1080
            sca
            clc
            disp('wrong resolution')
        end
        cd('C:\Users\PsychEEGlab\Desktop\EXPERIMENTS\Matt\Repetition supression\Data')
    end
    %% makes the condition table
    % This makes a condition table where the orientations are cued either
    % validly or invalidly with 75% accuracy in different blocks.
    % Furthermore, it makes 20% of the trials contain targets that the
    % subject has to respond to. The counter-balancing order is between
    % blocks with repetitions or alternations.
    if CountBalance ==1
        start=1;
    else
        start=2;
    end
    nblocks=20;
    ConditionTable=zeros(nblocks*135,4);
    count=0;
    nreps=3;
    for i =1:nblocks
        Temp1=Conditions(nreps,1,9);
        Temp2=Conditions(nreps,1,9);
        Temp3=Conditions(nreps,1,9);
        Temp4=Conditions(nreps,1,9);
        Temp4(:,2)=2;
        Temp=zeros(length(Temp4)*5,5);
        Temp(1:length(Temp4)*4,1:3)=[Temp1; Temp2 ;Temp3; Temp4];
        Temp(1:length(Temp4)*4,5)=1;
        Temp(length(Temp4)*4+1:end,5)=2;
        Temp5=Conditions(nreps,1,9);
        Temp(length(Temp4)*4+1:end,1:3)=Temp5;
        
        
        order=randperm(length(Temp));
        Temp=Temp(order,:);
        ConditionTable(count+1:count+length(Temp),1:5)=Temp;
        ConditionTable(count+1:count+length(Temp),4)=mod(start+i,2)+1;
        count=count+length(Temp);
    end
    ConditionTable(:,1)=1:length(ConditionTable);
    ResponseStore=zeros(length(ConditionTable),18);
    BlockLength=length(Temp);
    %% makes the orientations of T1
    T1ExpectedAlt=Conditions(nblocks*4.5,1,9);
    T1UnexpectedAlt=Conditions(nblocks*1.5,1,9);
    TargetAlt=Conditions(nblocks*1.5,1,9);
    
    ExpectedAcount=1;
    UnexpectedAcount=1;
    TargetAltcount=1;
    
    Start=GetSecs;
    for trial =1:length(ConditionTable)
        % expected repetition blocks
        if ConditionTable(trial,5) ==1
            if ConditionTable(trial,4)==1 % T1 == T2
                if ConditionTable(trial,2) ==1
                    ConditionTable(trial,6)=ConditionTable(trial,3);
                elseif ConditionTable(trial,2) ==2 % T2 =~ T2
                    while ConditionTable(trial,3)==T1UnexpectedAlt(UnexpectedAcount,3);
                        T1UnexpectedAlt(end+1,:)=T1UnexpectedAlt(UnexpectedAcount,:);
                        UnexpectedAcount=UnexpectedAcount+1;
                        if GetSecs-Start > 2
                            clc
                            disp('Try again : condition table broken')
                            sca
                            break
                            
                        end
                    end
                    
                    ConditionTable(trial,6)=T1UnexpectedAlt(UnexpectedAcount,3);
                    UnexpectedAcount=UnexpectedAcount+1;
                end
                % expected alternation block
            elseif ConditionTable(trial,4)==2 %T1 ~= T2
                if ConditionTable(trial,2) ==1
                    while ConditionTable(trial,3)==T1ExpectedAlt(ExpectedAcount,3);
                        T1ExpectedAlt(end+1,:)=T1ExpectedAlt(ExpectedAcount,:);
                        ExpectedAcount=ExpectedAcount+1;
                        if GetSecs-Start > 2
                            clc
                            disp('Try again : condition table broken')
                            sca
                            break
                        end
                    end
                    
                    ConditionTable(trial,6)=T1ExpectedAlt(ExpectedAcount,3);
                    ExpectedAcount=ExpectedAcount+1;
                    
                    
                else
                    ConditionTable(trial,6)=ConditionTable(trial,3);
                end
            end
            
        else
            if ConditionTable(trial,4) ==1
                ConditionTable(trial,6) =ConditionTable(trial,3);
            else
                ConditionTable(trial,6) =TargetAlt(TargetAltcount,3);
                TargetAltcount=TargetAltcount+1;
                
            end
        end
    end
    csvwrite(dataFileNameBU,ConditionTable);
    %% sets up the response keys
    RShift = KbName('RightShift');
    LShift = KbName('LeftShift');
    
    %% Select specific text font, style and size:
    Screen('TextFont',w, 'Helvetica');
    Screen('TextSize',w, 24);
    %% experimental parameters
    totframes=Hz*1.1; % total number of frames (1.25 secs x refresh rate)
    Idim=340; % Size of the overall gratings (each is half this)
    sc = Idim/7; % sigma of the guassian function
    freq = 18/Idim; % spatial frequency (cycles per pixel)
    contrast = 50; % contrast of the Gabor
    aspectratio = 1.0;
    gabortex = CreateProceduralGabor(w, Idim, Idim, 0, [0.5 0.5 0.5 0.0]);
    %% fixation points parameters
    fixLum= [0 0 0];
    fixSize=6;
    fixType=1;
    %%
    TargetTrials=sum(ConditionTable(:,5)==2)+4;
    TimeStore=zeros(totframes,length(ConditionTable)); % where stimulus presentation times will be stored
    AdaptorTime=1:12; % When the adaptor will be presented
    TestTime=73:84;% When the test will be presented
    WaitTimes=.1+rand(1,ntrials)*.2; % the inter-trial interval
    FirstTarget=randperm(TargetTrials)<TargetTrials/2;
    GreenTarget=randperm(TargetTrials)<TargetTrials/2;
    blockCounter=1; % block counter for the feedback screen
    TargetTrial=1;
    
    %%  Gives the person some instructions
    Screen('DrawTexture', w, gabortex, [], [center(1)-Idim/2-200 center(2)-Idim/2-250 center(1)+Idim/2-200 center(2)+Idim/2-250],  90, [], [], [1 RedMod RedMod], [],...
        kPsychDontDoRotation, [0, freq, sc, contrast, aspectratio, 0, 0, 0]);
    Screen('DrawTexture', w, gabortex, [],[center(1)-Idim/2+200 center(2)-Idim/2-250 center(1)+Idim/2+200 center(2)+Idim/2-250],  90, [], [], [ GreenMod 1 GreenMod], [],...
        kPsychDontDoRotation, [0, freq, sc, contrast, aspectratio, 0, 0, 0]);
    DrawFormattedText(w, 'Left shift key for red gratings \n\n \n\n Right shift key for green gratings ', 'center', 'center', 1);
    Screen('Flip', w);
    KbWait;
    Screen('Flip', w);
    
    WaitSecs(.5);
    
    %% waits until subject is ready to start
    DrawFormattedText(w, 'Press any key when ready to go', 'center', 'center', 1);
    Screen('Flip', w);
    KbWait;
    Screen('Flip', w);
    WaitSecs(.5);
    %% the trial loop
    Lasttrial=GetSecs;
    
    for trial =1:length(ConditionTable)
        %% gets the orientaitons for the two gratings for each trial
        TestOrientation=ConditionTable(trial,3)*20-20;
        AdaptOrientation=ConditionTable(trial,6)*20-20;
        
        if ConditionTable(trial,5)==1
            modulateColor1=[1 1 1];
            modulateColor2=[1 1 1];
        else
            if GreenTarget(TargetTrial) ==0
                if FirstTarget(TargetTrial) ==1
                    modulateColor1=[1 RedMod RedMod];
                    modulateColor2=[1 1 1];
                else
                    modulateColor1=[1 1 1];
                    modulateColor2=[1 RedMod RedMod];
                end
            else
                if FirstTarget(TargetTrial) ==1
                    modulateColor1=[GreenMod 1 GreenMod];
                    modulateColor2=[1 1 1];
                else
                    modulateColor1=[1 1 1];
                    modulateColor2=[GreenMod 1 GreenMod];
                end
            end
            TargetTrial=TargetTrial+1;
        end
        TrialPhase=randperm(180,1);
        TestPhase=randperm(180,1);
        
        %% makes the triggers
        trigger_value=ConditionTable(trial,3);
        
        if ConditionTable(trial,2)==2
            trigger_value=trigger_value+10;
        end
        if ConditionTable(trial,4)==2
            trigger_value=trigger_value+50;
        end
        
        if ConditionTable(trial,5)==2
            trigger_value=trigger_value+100;
        end
        %% waits until the required intertrial wait time
        while GetSecs-Lasttrial <WaitTimes(trial)
            WaitSecs(.001); % stops the computer from getting overloaded
        end
        KbQueueCreate;
        keyIsDown=0;
        
        %% does a pre-flip for timing issues (200 ms or 500 ms in first trial )
        if trial ==1
            for i=1:30
                Screen('DrawDots', w, [center(1) center(2)],fixSize,fixLum,[],fixType);
                Screen('DrawingFinished', w);
                Screen('Flip', w);
            end
        else
            for i=1:12
                Screen('DrawDots', w, [center(1) center(2)],fixSize,fixLum,[],fixType);
                Screen('DrawingFinished', w);
                Screen('Flip', w);
            end
        end
        KbQueueStart;
        
        %% send the trigger
        if send_triggers
            io64(ioObj,parallel_port_address,trigger_value); % set the port called parallel_port_address to your desired trigger value
            wait(pulse_duration); % wait a while
            io64(ioObj,parallel_port_address,0); % set the port called parallel_port_address back to zero
        end
        %% does the stimulus presentation
        vbl=Screen('Flip', w);
        startTime = GetSecs;
        for frame =1:totframes
            if ismember(frame,AdaptorTime)
                Screen('DrawTexture', w, gabortex, [], [], 90+AdaptOrientation, [], [], modulateColor1, [],...
                    kPsychDontDoRotation, [TrialPhase, freq, sc, contrast, aspectratio, 0, 0, 0]);
            end
            if ismember(frame,TestTime)
                Screen('DrawTexture', w, gabortex, [], [], 90+TestOrientation, [], [], modulateColor2, [],...
                    kPsychDontDoRotation, [TestPhase, freq, sc, contrast, aspectratio, 0, 0, 0]);
            end

            Screen('DrawingFinished', w);
            vbl=Screen('Flip', w, vbl + (1-0.5)*ifi);
            TimeStore(frame,trial)=vbl;
        end
        %%
        Lasttrial=GetSecs; % end of trial timer
        TrialTime=Lasttrial-startTime;
        Screen('DrawDots', w, [center(1) center(2)],fixSize,fixLum,[],fixType);
        Screen('DrawingFinished', w);
        Screen('Flip', w, vbl + (1-0.5)*ifi);
        %% checks to see whether keys were pressed
        [keyIsDown, FirstPress]=KbQueueCheck;
        KbQueueRelease;
        KbQueueStop;
        %% gets respones
        if ConditionTable(trial,5) ==1
            if keyIsDown==0
                RT=NaN;
                correct=98; % correct rejection
            else
                correct=99; % false alarm
                RT=min((FirstPress(logical(FirstPress))-startTime)*1000);
            end
        else
            if keyIsDown==1
                if GreenTarget(TargetTrial-1) ==1
                    if (FirstPress(RShift))
                        correct=1; % correct response
                    elseif (FirstPress(LShift))
                        correct=0; % incorrect response
                    end
                else
                    if (FirstPress(RShift))
                        correct=0; % incorrect response
                    elseif (FirstPress(LShift))
                        correct=1; % correct response
                    end
                end
                if FirstTarget(TargetTrial-1) ==1
                    RT=min((FirstPress(logical(FirstPress))-startTime)*1000);
                else
                    RT=min((FirstPress(logical(FirstPress))-startTime)*1000)-500;
                end
            else
                correct=3; % missed
                RT=NaN;
            end
        end
        %% sends trigger about correct/incorrect response
        if correct ==1;
            Response_trigger=240;
        elseif correct==0;
            Response_trigger=241;
        elseif correct==99
            Response_trigger=242;
        else
            Response_trigger=243;
        end
        if send_triggers
            io64(ioObj,parallel_port_address,Response_trigger); % set the port called parallel_port_address to your desired trigger value
            wait(pulse_duration); % wait a while
            io64(ioObj,parallel_port_address,0); % set the port called parallel_port_address back to zero
        end
        FirstPress={};
        %% puts the responses into the data store
        ResponseStore(trial,:)=[SubNo,trial,ConditionTable(trial,2),...
            ConditionTable(trial,3), ConditionTable(trial,4), ...
            ConditionTable(trial,5), WaitTimes(trial) , GreenTarget(TargetTrial), ...
            FirstTarget(TargetTrial),TrialPhase,TestPhase,AdaptOrientation,TestOrientation,...
            trigger_value,Response_trigger, TrialTime, RT, correct]  ;
        %% writes the responses to a file
        dataFile = fopen(dataFileName, 'a');
        fprintf(dataFile, formatString,ResponseStore(trial,:));
        fclose(dataFile);
        %% displays the trial information to the command window
        formatSpec = 'Block : %d | Trial : %d | Trial time :  %4.3f | Trigger : %d \n';
        fprintf(formatSpec,blockCounter, trial, TrialTime,trigger_value)
        %% asks the subject to have a break every X many trials
        if mod(trial,BlockLength) ==0 && trial ~=length(ConditionTable)
            blocktrials=trial-BlockLength+1:trial;
            TT=ResponseStore(blocktrials,6)==2;
            ResponsesBlock=ResponseStore(blocktrials,18);
            BlockCorrect=mean(ResponsesBlock(TT)==1);
            MissedBlock=mean(ResponsesBlock(TT)==3);
            AvgRT=ResponseStore(blocktrials,17);
            AvgRT=nanmean(AvgRT(TT));
            FA=sum(ResponseStore(blocktrials,18)==99);
            %% gives the subject feedback about performance in the last block
            myText=strcat('Block done (', num2str(blockCounter),'/',num2str(nblocks), ')\n \n\n',...
                'Correct :',num2str(round(BlockCorrect*100)),'% \n \n',...
                'Missed :' ,num2str(round(MissedBlock*100)),'%  \n \n',...
                'False alarms :  ',num2str(FA), ' \n\n' ,...
                'Median reaction time :  ',num2str(round(AvgRT)), ' ms \n\n \n\n' ,...
                'Press any key when ready to go');
            DrawFormattedText(w, myText,  'center', 'center', 1);
            Screen('Flip', w); % flips the feedback to the screen
            %% displays the trial information to the command window
            formatSpec = 'Block over - Correct : %d | Missed :  %d | False alarms : %d \n';
            fprintf(formatSpec, round(BlockCorrect*100),round(MissedBlock*100), round(MissedBlock*100),FA)
            %%
            WaitSecs(1);
            blockCounter=blockCounter+1; % increments the block counter
            %% waits for the subject to initiate the next block
            KbWait;
            Screen('Flip', w); % flips a blank screen
            WaitSecs(.2);
            %%  Gives the person some instructions
            Screen('DrawTexture', w, gabortex, [], [center(1)-Idim/2-200 center(2)-Idim/2-250 center(1)+Idim/2-200 center(2)+Idim/2-250],  90, [], [], [1 RedMod RedMod], [],...
                kPsychDontDoRotation, [0, freq, sc, contrast, aspectratio, 0, 0, 0]);
            Screen('DrawTexture', w, gabortex, [],[center(1)-Idim/2+200 center(2)-Idim/2-250 center(1)+Idim/2+200 center(2)+Idim/2-250],  90, [], [], [ GreenMod 1 GreenMod], [],...
                kPsychDontDoRotation, [0, freq, sc, contrast, aspectratio, 0, 0, 0]);
            DrawFormattedText(w, 'Left shift key for red gratings \n\n \n\n Right shift key for green gratings ', 'center', 'center', 1);
            Screen('Flip', w);
            KbWait;
            Screen('Flip', w);
            WaitSecs(1);
        end
    end
catch
    ple
    sca
end
sca

plot(diff(TimeStore));