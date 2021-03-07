%% Experiment Script
clear all;
clc;
close all;


%%%%%%%%%%%%%%%%%%%%%%%%% ID %%%%%%%%%%%%%%%%%%%%%%%%%%%%

ID = '';
while (length(ID) == 0) || (ID <= -1);
    ID = input('Enter participant ID: ');
    if length(ID)==0       

            fprintf('No ID was provided. \n');
    
    elseif ID <= -1
            fprintf ('Invalid ID. \n');
        
    else  
        fprintf('Participant ID = %d\n \n', ID);
    end 
end 

%%%%%%%%%%%%%%%%%%%%% Number of Trials %%%%%%%%%%%%%%%%%%%%%%%%%%%

n = '';
while (isempty(n)) || (n <= -1) || mod(n,2)==1
    n = input('Please enter a positive, even number of trials: ');
    if isempty(n)
        fprintf('No number was provided! \n');
        
    elseif n <= -1
        fprintf('Invalid number was provided! \n');
   
    elseif isnan(n)
        fprintf('Invalid input was provided! Please enter a numerical value : \n');
    
    elseif mod(n,2)==1
        fprintf('An odd number was provided! \n');
    end
end 

fprintf('Number of trials = %d\n \n', n);

%%%%%%%%%%%%%%%%%%%%% Standardize Randomization %%%%%%%%%%%%%%%%%%%%%%%%%%


rng(ID); %Seed

conditions = [ones(1, n/2) zeros(1, n/2)];

r = randperm(n);

conditions = conditions(r); %conditional array has now been randomized, but the number of 1s and 0s remaisn the same as before
    
%0 is noise; 1 is signal.

%%%%%%%%%%%%%%%%%%%%%%%%% Sound Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sr = 44100; %the sample rate in Hz
gatedur = .01; %ramp for smoothing sound
gate = cos(linspace(pi, 2*pi, sr*gatedur));
gate = gate+1;
gate = gate/2;
offsetgate = fliplr(gate);
d = 0.05;
t = linspace(0, d, sr*d);

frequencies = [];
signal_tones = [];
noise_tones = [];

%Singal setup
for signal = 1:4
        sequence = [];
        f = 100.*randn(1) + 1100; 
        frequencies = [frequencies, f];
        tone = sin(2*pi*f*t);
        sustain = ones(1, (length(tone)-2*length(gate)));
        envelope = [gate, sustain, offsetgate];
        tone = envelope .* tone;
        silence = zeros(1, sr*d); 
        signal_tones = [signal_tones, tone, silence];
end 

%Noise setup
for noise = 1:4
        sequence = [];
        f = 10.*randn(1) + 1000;
        frequencies = [frequencies, f];
        tone = sin(2*pi*f*t);
        sustain = ones(1, (length(tone)-2*length(gate)));
        envelope = [gate, sustain, offsetgate];
        tone = envelope .* tone;
        silence = zeros(1, sr*d);
        noise_tones = [noise_tones, tone, silence];
end 

% The specific frequencies is assigned in the trial interface below to the
% appropreiate conditions

%%%%%%%%%%%%%%%%%%%%%%%%% Trail Interface %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Section not finished - still a few bugs

data = [];
hits = 0;
FAs = 0;
% while ((data == 0) & (condition == 1)) | ((data == 1) & (condition == 0)) | (length(data) == 0) | (data <= -1) | (data >= 2);
    for i = 1:n
        fprintf('Trial %d of %d \n', i, n);
        condition = conditions(i);

        %Enable the following line to see the correct answer each time;
        %disable this line when submitting assignment
        %fprintf('Trial %d condition = %d\n', i , condition);
        
        %Enable the following lines to make the program play sound
        if condition == 1 %If condition is signal
            soundsc([signal_tones], sr); %Then play signal melody
        
        elseif condition == 0 %If condition is noise
            soundsc([noise_tones], sr); %Then play noise melody
        end 
        
        while (length(data) == 0) | (data <= -1) | (data >= 2) | (data == 0) | (data == 1)
        data = input('Was this a signal (1) or noise (0)?: ');

            if (length(data) == 0) | (data <= -1) | (data >= 2)
                fprintf('\tInvalid input! Please provide a valid input. \n \n');
                continue
            elseif (data == 0) | (data == 1)
%                 continue
                break
            end
        end
            
%             continue 
            
%             break
        %Under these conditions, the program will restart from trial 1
        %instead of remaining within the same trial with "break", or will
        %advance to the next trial without "break"
        
        %I have tried to put a while loop within the condition but it
        %resulted in an infinite loop - Josh

        if ((data == 0) & (condition == 1)) 
            fprintf('\tIncorrect! \n \n')
            
            continue
            
        elseif ((data == 1) & (condition == 0))
            fprintf('\tIncorrect! \n \n')
            FAs = FAs + 1;
            continue

        elseif ((data == 1) & (condition == 1)) 
            fprintf('\tCorrect! \n \n')
            hits = hits + 1;
            continue
            
        elseif ((data == 0) & (condition == 0))
            fprintf('\tCorrect! \n \n')
            
            continue 
        end
    end
% end
% I have disabled the while-loop to force the trials to move on - Josh

% display(data)
% display(hits)
% display(FAs)

%%%%%%%%%%%%%%%%%%%%% Calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if rem(n, 2) == 1
    n = n + 1
end

fprintf('Data of Participant %d: \n', ID)

%Prob of hit
pHit = hits / (n/2);
if pHit == 1;
    pHit = 0.99;
elseif pHit == 0;
    pHit = 0.01;
end 

fprintf('\tProbablity of hit: %0.2f. \n', pHit);

%Prob of false alarm
pFA = FAs/ (n/2);
if pFA == 1;
    pFA = 0.99;
elseif pFA == 0;
    pFA = 0.01;
end 

fprintf('\tProbablity of false alarm: %0.2f. \n', pFA);

%Sensitivity
zHit = -sqrt(2)*erfcinv(2*pHit);
zFA = -sqrt(2)*erfcinv(2*pFA);
d = zHit-zFA;

fprintf('\tSensitivity: %2.1f \n', d);

%Bias index
B = exp(((-sqrt(2)*erfcinv(2*pHit))^2 - (-sqrt(2)*erfcinv(2*pFA))^2)/2);

fprintf('\tBias index: %2.1f \n', B);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Saving the File %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

participant_data = sprintf('Participant %d data.mat', ID);
save(participant_data, 'ID', 'n', 'r', 'conditions', 'hits', 'FAs', 'pHit', 'pFA', 'd', 'B')

%%%%%%%%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%listing all .mat files in this name; '*' being participant No.
% X = struct2table(dir('Participant * data.mat'))

% load('Participant 7 data.mat');







