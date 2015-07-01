try
%% Setup
stimsize = 170; % size in pix
transition_duration = 1.5; % duration in seconds
transition_frames = round(transition_duration*60); % frames (adjust for screen)

rng('shuffle');

swap=@(varargin)varargin{nargin:-1:1}; % this function is used later to swap the variables

% this is the time the transitions are going to take on the criterion
% trials
mean_dur_contr = 2;
mean_dur_asc = 3;

while ~exist('participant_id', 'var') || isempty(participant_id)
    participant_id = inputdlg('Please enter a participant id:');
end

secs = zeros(100000, 1); sacc_frames = cell(20000, 1); keypresses = zeros(100000, 3); %pre-allocate for speed


KbName('UnifyKeynames');
arrowkeys = [KbName('LeftArrow'), KbName('UpArrow'), KbName('RightArrow')];



screen.no = 2;
handle = Screen('OpenWindow', 2, [255 215 0]);
screen.dimensions=Screen('Rect', handle);
screen.xcen = screen.dimensions(3)/2;
screen.ycen = screen.dimensions(4)/2;
Screen('BlendFunction', handle, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
Screen('DrawText', handle, 'Preparing the stimuli. Please wait.');
Screen('Flip', handle);
HideCursor;

%% Stimulus positions
d = stimsize; % this will have to come from calibration later

circle_rects = cat(2,([screen.xcen-stimsize/2+d, screen.ycen-stimsize/2, screen.xcen+stimsize/2+d, screen.ycen+stimsize/2] + [-6, -6, +6, +6])',...
                        ([screen.xcen-stimsize/2-d, screen.ycen-stimsize/2, screen.xcen+stimsize/2-d, screen.ycen+stimsize/2] + [-6, -6, +6, +6])');

target_rect_right = circle_rects(:, 1)' + [(sqrt((stimsize^2)/2)/4) (sqrt((stimsize^2)/2)/4) -(sqrt((stimsize^2)/2)/4) -(sqrt((stimsize^2)/2)/4)];
target_rect_left = circle_rects(:, 2)'  + [(sqrt((stimsize^2)/2)/4) (sqrt((stimsize^2)/2)/4) -(sqrt((stimsize^2)/2)/4) -(sqrt((stimsize^2)/2)/4)];

%% making of textures
% These are the images (should be randomised between trials)
pic1 = ['D:\Dropbox\Documents\University\Matlab Code\Rivalry\Simulation\Pictures\green' num2str(1) '.JPG'];
pic2 = ['D:\Dropbox\Documents\University\Matlab Code\Rivalry\Simulation\Pictures\red' num2str(1) '.JPG'];

pic1 = imread(pic1); pic1 = imresize(pic1, stimsize/size(pic1, 1));
pic2 = imread(pic2); pic2 = imresize(pic2, stimsize/size(pic2, 1));

pic_texture(1) = Screen('MakeTexture', handle, pic1);
pic_texture(2) = Screen('MakeTexture', handle, pic2);

n_tex = 5;

for i = 1:n_tex

eval(['load gaussians', num2str(i), '.mat']); % there are 10 different gauss files, so we load them one after the other
for ii = 1:transition_frames
    tex = round(linspace(1, 180, transition_frames));
    tex = tex(ii);
    curr_gauss = imresize(gaussians(tex).red, stimsize/500);
    pic_textures(1, ii, i) = Screen('MakeTexture', handle, cat(3, pic1, curr_gauss));
    pic_textures(2, ii, i) = Screen('MakeTexture', handle, cat(3, pic2, curr_gauss));
end
end
clear gaussians curr_gauss tex 

%% Find the rivalry spot
while 1
    circle_rects = cat(2,([screen.xcen-stimsize/2+d, screen.ycen-stimsize/2, screen.xcen+stimsize/2+d, screen.ycen+stimsize/2] + [-6, -6, +6, +6])',...
                            ([screen.xcen-stimsize/2-d, screen.ycen-stimsize/2, screen.xcen+stimsize/2-d, screen.ycen+stimsize/2] + [-6, -6, +6, +6])');

    target_rect_right = circle_rects(:, 1)' + [(sqrt((stimsize^2)/2)/4) (sqrt((stimsize^2)/2)/4) -(sqrt((stimsize^2)/2)/4) -(sqrt((stimsize^2)/2)/4)];
    target_rect_left = circle_rects(:, 2)'  + [(sqrt((stimsize^2)/2)/4) (sqrt((stimsize^2)/2)/4) -(sqrt((stimsize^2)/2)/4) -(sqrt((stimsize^2)/2)/4)];
    
    
    Screen('FrameOval', handle, [0 0 0], circle_rects, 6);
    Screen('Flip', handle);
    WaitSecs(0.2);
    [~, keyCode, ~] = KbWait;
    
    if keyCode(KbName('UpArrow'))
        d = d+4;
    elseif keyCode(KbName('DownArrow'))
        d = d-4;
    elseif keyCode(KbName('Return'))
        d = d+stimsize/2;
    elseif keyCode(KbName('Escape'))
        break
    end
    
end





%% Demonstrate the stimuli
% red image
Screen('DrawTexture', handle, pic_texture(2), [], target_rect_right);
Screen('DrawTexture', handle, pic_texture(2), [], target_rect_left);
Screen('FrameOval', handle, [0 0 0], circle_rects, 6);
Screen('Flip', handle);
WaitSecs(1);
KbWait;
% green image
Screen('DrawTexture', handle, pic_texture(1), [], target_rect_right);
Screen('DrawTexture', handle, pic_texture(1), [], target_rect_left);
Screen('FrameOval', handle, [0 0 0], circle_rects, 6);
Screen('Flip', handle);
WaitSecs(1);
KbWait;
% mixed image
Screen('DrawTexture', handle, pic_texture(2), [], target_rect_right);
Screen('DrawTexture', handle, pic_texture(2), [], target_rect_left);
Screen('DrawTexture', handle, pic_textures(1, transition_frames/2, 1), [], target_rect_right);
Screen('DrawTexture', handle, pic_textures(1, transition_frames/2, 1), [], target_rect_left);
Screen('FrameOval', handle, [0 0 0], circle_rects, 6);
Screen('Flip', handle);
WaitSecs(1);
KbWait;

Screen('FrameOval', handle, [0 0 0], circle_rects, 6);
Screen('Flip', handle);
WaitSecs(1);
KbWait;

%% Pracice trial: Slow smooth transitions
for trialcount = 1:2
%% Prepare Schedule
j = 1;
load('event_durations.mat');
total_time = 0;
if trialcount == 1 || trialcount == 2
% values for controls
    mu_switches = 17.416;
    sd_switches = 5.01;
    mu_reversions = 1.736;
    sd_reversions = 1.19;
    DomDistr = adjDomC;
    MixDistr = adjMixC;
    meanMix = meanMixC;
elseif trialcount == 3 || trialcount == 4
    mu_switches = 17.416;
    sd_switches = 5.01;
    mu_reversions = 1.736;
    sd_reversions = 1.19;
    DomDistr = adjDomA;
    MixDistr = adjMixA;
    meanMix = meanMixA;
end


while total_time < 55 || total_time > 60
    clear transitions;
    no_of_switches = round(normrnd(mu_switches, sd_switches));
    while no_of_switches > mu_switches+sd_switches || no_of_switches < mu_switches-sd_switches
        no_of_switches = round(normrnd(mu_switches, sd_switches));
    end
    no_of_reversions = round(normrnd(mu_reversions, sd_reversions));
    while no_of_reversions > mu_reversions+sd_reversions || no_of_reversions < mu_reversions-sd_reversions
        no_of_reversions = round(normrnd(mu_reversions, sd_reversions));
    end
    
    no_of_doms = no_of_switches + no_of_reversions +1;

    dom_durations = DomDistr(randi(size(DomDistr, 1), 1, no_of_doms));
    mix_durations = MixDistr(randi(size(MixDistr, 1), 1, no_of_reversions+no_of_switches));
    transitions.type(1, 1:no_of_switches) = 1;   transitions.type(1, no_of_switches+1:no_of_switches+no_of_reversions) = 2;
    for i = 1:no_of_switches+no_of_reversions
        transitions.onset(1, i) = sum(mix_durations(1:i-1))+sum(dom_durations(1:i)); % this is the onset of this transition
        transitions.end(1, i) = sum(mix_durations(1:i))+sum(dom_durations(1:i)); % this is the end of this transition
    end

    warmup = rand*(60-transitions.end(1, no_of_switches+no_of_reversions));
    transitions.onset = transitions.onset + warmup;
    transitions.end = transitions.end + warmup;
    total_time = transitions.end(1, no_of_switches+no_of_reversions);
end
    mix_frames = round(mix_durations*60);
    transitions.type = transitions.type(randperm(no_of_switches+no_of_reversions));
    
    if trialcount == 1
    % just for the practice:
    clear transitions;
    transitions.type = [1 1 1 2 1 2];
    transitions.onset = [5 12 19 25 28 31];
    transitions.end = [8 16 23 26 30 35];
    mix_frames = (transitions.end - transitions.onset) * 60;
    warmup = 1;
    end
    
%% Fixed duration of transitions
Screen('FrameOval', handle, [0 0 0], circle_rects, 6);
Screen('Flip', handle);
KbWait;
t0 = GetSecs;

if round(rand)
    % first image will be red
    tex = randi(n_tex);
    frames = warmup*60;
    for i = round(linspace(transition_frames/2, transition_frames, frames))
    Screen('DrawTexture', handle, pic_texture(2), [], target_rect_right);
    Screen('DrawTexture', handle, pic_texture(2), [], target_rect_left);
    Screen('DrawTexture', handle, pic_textures(1, i, tex), [], target_rect_right);
    Screen('DrawTexture', handle, pic_textures(1, i, tex), [], target_rect_left);
    Screen('FrameOval', handle, [0 0 0], circle_rects, 6);
    Screen('Flip', handle);
    j=j+1; [~, secs(j, 1), keyCode, ~] = KbCheck; keypresses(j, :) = keyCode(arrowkeys); % sacc_frames{j, 1} = GetSaccadometerNewFrame(pstr);
    j=j+1; [~, secs(j, 1), keyCode, ~] = KbCheck; keypresses(j, :) = keyCode(arrowkeys); % sacc_frames{j, 1} = GetSaccadometerNewFrame(pstr);
    end
    last_image = 1;
    next_image = 2;
    
else
    % first image will be green
    tex = randi(n_tex);
    frames = warmup*60;
    for i = round(linspace(transition_frames/2, transition_frames, frames))
    Screen('DrawTexture', handle, pic_texture(1), [], target_rect_right);
    Screen('DrawTexture', handle, pic_texture(1), [], target_rect_left);
    Screen('DrawTexture', handle, pic_textures(2, i, tex), [], target_rect_right);
    Screen('DrawTexture', handle, pic_textures(2, i, tex), [], target_rect_left);
    Screen('FrameOval', handle, [0 0 0], circle_rects, 6);
    Screen('Flip', handle);
    j=j+1; [~, secs(j, 1), keyCode, ~] = KbCheck; keypresses(j, :) = keyCode(arrowkeys); % sacc_frames{j, 1} = GetSaccadometerNewFrame(pstr);
    j=j+1; [~, secs(j, 1), keyCode, ~] = KbCheck; keypresses(j, :) = keyCode(arrowkeys); % sacc_frames{j, 1} = GetSaccadometerNewFrame(pstr);
    end
    last_image = 2;
    next_image = 1;

end % we start on an ambiguous image


for ii = 1:size(transitions.type,2)
    
    while GetSecs < t0+transitions.onset(ii)
    j=j+1;
    [~, secs(j, 1), keyCode, ~] = KbCheck;
    %        sacc_frames{j, 1} = GetSaccadometerNewFrame(pstr);
    keypresses(j, :) = keyCode(arrowkeys); % 38-40: arrows, 68: d, 87: w
    end % we wait out the dominance period
    
    tex = randi(n_tex);

    if transitions.type(ii) == 1 % switch
        
        for i = round(linspace(1, transition_frames, mix_frames(ii)))
        Screen('DrawTexture', handle, pic_texture(last_image), [], target_rect_right);
        Screen('DrawTexture', handle, pic_texture(last_image), [], target_rect_left);
        Screen('DrawTexture', handle, pic_textures(next_image, i, tex), [], target_rect_right);
        Screen('DrawTexture', handle, pic_textures(next_image, i, tex), [], target_rect_left);
        Screen('FrameOval', handle, [0 0 0], circle_rects, 6);
        Screen('Flip', handle);
        j=j+1; [~, secs(j, 1), keyCode, ~] = KbCheck; keypresses(j, :) = keyCode(arrowkeys); % sacc_frames{j, 1} = GetSaccadometerNewFrame(pstr);
        j=j+1; [~, secs(j, 1), keyCode, ~] = KbCheck; keypresses(j, :) = keyCode(arrowkeys); % sacc_frames{j, 1} = GetSaccadometerNewFrame(pstr);
        end
        
        [next_image, last_image] = swap(next_image, last_image);
        
    elseif transitions.type(ii) == 2 % reversion
        
        for i = [round(linspace(1, transition_frames/2, mix_frames(ii)/2)), round(linspace(transition_frames/2, 1, mix_frames(ii)/2))]
        Screen('DrawTexture', handle, pic_texture(last_image), [], target_rect_right);
        Screen('DrawTexture', handle, pic_texture(last_image), [], target_rect_left);
        Screen('DrawTexture', handle, pic_textures(next_image, i, tex), [], target_rect_right);
        Screen('DrawTexture', handle, pic_textures(next_image, i, tex), [], target_rect_left);
        Screen('FrameOval', handle, [0 0 0], circle_rects, 6);
        Screen('Flip', handle);
        j=j+1; [~, secs(j, 1), keyCode, ~] = KbCheck; keypresses(j, :) = keyCode(arrowkeys); % sacc_frames{j, 1} = GetSaccadometerNewFrame(pstr);
        j=j+1; [~, secs(j, 1), keyCode, ~] = KbCheck; keypresses(j, :) = keyCode(arrowkeys); % sacc_frames{j, 1} = GetSaccadometerNewFrame(pstr);
        end
        
    end
    
    

    
end
    
while GetSecs < t0+36
j=j+1;
[~, secs(j, 1), keyCode, ~] = KbCheck;
%        sacc_frames{j, 1} = GetSaccadometerNewFrame(pstr);
keypresses(j, :) = keyCode(arrowkeys); % 38-40: arrows, 68: d, 87: w
end % we wait until 60s is reached
Screen('Flip', handle);

%% Reset variables
secs = zeros(100000, 1); sacc_frames = cell(20000, 1); keypresses = zeros(100000, 3);
j = 1;

%% make new textures
pic1 = ['D:\Dropbox\Documents\University\Matlab Code\Rivalry\Simulation\Pictures\green' num2str(trialcount+1) '.JPG'];
pic2 = ['D:\Dropbox\Documents\University\Matlab Code\Rivalry\Simulation\Pictures\red' num2str(trialcount+1) '.JPG'];

pic1 = imread(pic1); pic1 = imresize(pic1, stimsize/size(pic1, 1));
pic2 = imread(pic2); pic2 = imresize(pic2, stimsize/size(pic2, 1));

pic_texture(1) = Screen('MakeTexture', handle, pic1);
pic_texture(2) = Screen('MakeTexture', handle, pic2);

n_tex = 5;
for i = 1:n_tex

eval(['load gaussians', num2str(i), '.mat']); % there are 10 different gauss files, so we load them one after the other
for ii = 1:transition_frames
    tex = round(linspace(1, 180, transition_frames));
    tex = tex(ii);
    curr_gauss = imresize(gaussians(tex).red, stimsize/500);
    pic_textures(1, ii, i) = Screen('MakeTexture', handle, cat(3, pic1, curr_gauss));
    pic_textures(2, ii, i) = Screen('MakeTexture', handle, cat(3, pic2, curr_gauss));
end
end
clear gaussians curr_gauss tex

end



%% Trial loop: criterion trials
for trialcount = 1:4
j = 1;

%% Prepare Schedule
load('event_durations.mat');
total_time = 0;
clear transitions;
if trialcount == 1 || trialcount == 2
% values for controls
    mu_switches = 17.416;
    sd_switches = 5.01;
    mu_reversions = 1.736;
    sd_reversions = 1.19;
    DomDistr = adjDomC;
    MixDistr = adjMixC;
    meanMix = meanMixC;
    total_time = 0;
elseif trialcount == 3 || trialcount == 4
    mu_switches = 17.416;
    sd_switches = 5.01;
    mu_reversions = 1.736;
    sd_reversions = 1.19;
    DomDistr = adjDomA;
    MixDistr = adjMixA;
    meanMix = meanMixA;
    total_time = 0;
end


while total_time < 55 || total_time > 59
    no_of_switches = round(normrnd(mu_switches, sd_switches));
    while no_of_switches > mu_switches+sd_switches || no_of_switches < mu_switches-sd_switches
        no_of_switches = round(normrnd(mu_switches, sd_switches));
    end
    no_of_reversions = round(normrnd(mu_reversions, sd_reversions));
    while no_of_reversions > mu_reversions+sd_reversions || no_of_reversions < mu_reversions-sd_reversions
        no_of_reversions = round(normrnd(mu_reversions, sd_reversions));
    end
    
    no_of_doms = no_of_switches + no_of_reversions +1;

    dom_durations = DomDistr(randi(size(DomDistr, 1), 1, no_of_doms));
    mix_durations = ones(1, no_of_reversions+no_of_switches)*meanMix;
    total_time = sum(dom_durations) + sum(mix_durations);
    
    transitions.type(1, 1:no_of_switches) = 1;   transitions.type(1, no_of_switches+1:no_of_switches+no_of_reversions) = 2;
    for i = 1:no_of_switches+no_of_reversions
        transitions.onset(1, i) = sum(mix_durations(1:i-1))+sum(dom_durations(1:i)); % this is the onset of this transition
        transitions.end(1, i) = sum(mix_durations(1:i))+sum(dom_durations(1:i)); % this is the end of this transition
    end
    warmup = rand*(60-total_time);
    transitions.onset = transitions.onset + warmup;
    transitions.end = transitions.end + warmup;

end
    mix_frames = round(mix_durations*60);
    transitions.type = transitions.type(randperm(no_of_switches+no_of_reversions));

%% Fixed duration of transitions
Screen('FrameOval', handle, [0 0 0], circle_rects, 6);
Screen('Flip', handle);
KbWait;
t0 = GetSecs;

if round(rand)
    % first image will be red
    tex = randi(n_tex);
    frames = warmup*60;
    for i = round(linspace(transition_frames/2, transition_frames, frames))
    Screen('DrawTexture', handle, pic_texture(2), [], target_rect_right);
    Screen('DrawTexture', handle, pic_texture(2), [], target_rect_left);
    Screen('DrawTexture', handle, pic_textures(1, i, tex), [], target_rect_right);
    Screen('DrawTexture', handle, pic_textures(1, i, tex), [], target_rect_left);
    Screen('FrameOval', handle, [0 0 0], circle_rects, 6);
    Screen('Flip', handle);
    j=j+1; [~, secs(j, 1), keyCode, ~] = KbCheck; keypresses(j, :) = keyCode(arrowkeys); % sacc_frames{j, 1} = GetSaccadometerNewFrame(pstr);
    j=j+1; [~, secs(j, 1), keyCode, ~] = KbCheck; keypresses(j, :) = keyCode(arrowkeys); % sacc_frames{j, 1} = GetSaccadometerNewFrame(pstr);
    end
    first_dom = 'red';
    last_image = 1;
    next_image = 2;
    
else
    % first image will be green
    tex = randi(n_tex);
    frames = warmup*60;
    for i = round(linspace(transition_frames/2, transition_frames, frames))
    Screen('DrawTexture', handle, pic_texture(1), [], target_rect_right);
    Screen('DrawTexture', handle, pic_texture(1), [], target_rect_left);
    Screen('DrawTexture', handle, pic_textures(2, i, tex), [], target_rect_right);
    Screen('DrawTexture', handle, pic_textures(2, i, tex), [], target_rect_left);
    Screen('FrameOval', handle, [0 0 0], circle_rects, 6);
    Screen('Flip', handle);
    j=j+1; [~, secs(j, 1), keyCode, ~] = KbCheck; keypresses(j, :) = keyCode(arrowkeys); % sacc_frames{j, 1} = GetSaccadometerNewFrame(pstr);
    j=j+1; [~, secs(j, 1), keyCode, ~] = KbCheck; keypresses(j, :) = keyCode(arrowkeys); % sacc_frames{j, 1} = GetSaccadometerNewFrame(pstr);
    end
    first_dom = 'green';
    last_image = 2;
    next_image = 1;

end % we start on an ambiguous image


for ii = 1:no_of_switches+no_of_reversions
    
    while GetSecs < t0+transitions.onset(ii)
    j=j+1;
    [~, secs(j, 1), keyCode, ~] = KbCheck;
    %        sacc_frames{j, 1} = GetSaccadometerNewFrame(pstr);
    keypresses(j, :) = keyCode(arrowkeys); % 38-40: arrows, 68: d, 87: w
    end % we wait out the dominance period
    
    tex = randi(n_tex);

    if transitions.type(ii) == 1 % switch
        
        for i = round(linspace(1, transition_frames, mix_frames(ii)))
        Screen('DrawTexture', handle, pic_texture(last_image), [], target_rect_right);
        Screen('DrawTexture', handle, pic_texture(last_image), [], target_rect_left);
        Screen('DrawTexture', handle, pic_textures(next_image, i, tex), [], target_rect_right);
        Screen('DrawTexture', handle, pic_textures(next_image, i, tex), [], target_rect_left);
        Screen('FrameOval', handle, [0 0 0], circle_rects, 6);
        Screen('Flip', handle);
        j=j+1; [~, secs(j, 1), keyCode, ~] = KbCheck; keypresses(j, :) = keyCode(arrowkeys); % sacc_frames{j, 1} = GetSaccadometerNewFrame(pstr);
        j=j+1; [~, secs(j, 1), keyCode, ~] = KbCheck; keypresses(j, :) = keyCode(arrowkeys); % sacc_frames{j, 1} = GetSaccadometerNewFrame(pstr);
        end
        
        [next_image, last_image] = swap(next_image, last_image);
        
    elseif transitions.type(ii) == 2 % reversion
        
        for i = [round(linspace(1, transition_frames/2, mix_frames(ii)/2)), round(linspace(transition_frames/2, 1, mix_frames(ii)/2))]
        Screen('DrawTexture', handle, pic_texture(last_image), [], target_rect_right);
        Screen('DrawTexture', handle, pic_texture(last_image), [], target_rect_left);
        Screen('DrawTexture', handle, pic_textures(next_image, i, tex), [], target_rect_right);
        Screen('DrawTexture', handle, pic_textures(next_image, i, tex), [], target_rect_left);
        Screen('FrameOval', handle, [0 0 0], circle_rects, 6);
        Screen('Flip', handle);
        j=j+1; [~, secs(j, 1), keyCode, ~] = KbCheck; keypresses(j, :) = keyCode(arrowkeys); % sacc_frames{j, 1} = GetSaccadometerNewFrame(pstr);
        j=j+1; [~, secs(j, 1), keyCode, ~] = KbCheck; keypresses(j, :) = keyCode(arrowkeys); % sacc_frames{j, 1} = GetSaccadometerNewFrame(pstr);
        end
        
    end
    
    

    
end
    
while GetSecs < t0+60
j=j+1;
[~, secs(j, 1), keyCode, ~] = KbCheck;
%        sacc_frames{j, 1} = GetSaccadometerNewFrame(pstr);
keypresses(j, :) = keyCode(arrowkeys); % 38-40: arrows, 68: d, 87: w
end % we wait until 60s is reached
Screen('Flip', handle);
%% Save data and reset variables
filename = strcat('D:\Dropbox\Documents\University\Matlab Code\Rivalry\Simulation\Results\', num2str(date), '_', participant_id, '_simulation_CRITERION_trace_', num2str(trialcount), '.mat');
if iscell(filename)
    filename = filename{1};
end
secs = secs(1:j, 1); keypresses = keypresses(1:j, :); % sacc_frames = sacc_frames(1:j, 1);
secs = secs - t0;
save(filename, 'secs', 'keypresses', 'sacc_frames', 'transitions', 'mix_durations', 'dom_durations', 'first_dom');
secs = zeros(100000, 1); sacc_frames = cell(20000, 1); keypresses = zeros(100000, 3);
j = 1;

%% make new textures
pic1 = ['D:\Dropbox\Documents\University\Matlab Code\Rivalry\Simulation\Pictures\green' num2str(mod(trialcount+2, 6)+1) '.JPG'];
pic2 = ['D:\Dropbox\Documents\University\Matlab Code\Rivalry\Simulation\Pictures\red' num2str(mod(trialcount+2, 6)+1) '.JPG'];

pic1 = imread(pic1); pic1 = imresize(pic1, stimsize/size(pic1, 1));
pic2 = imread(pic2); pic2 = imresize(pic2, stimsize/size(pic2, 1));

pic_texture(1) = Screen('MakeTexture', handle, pic1);
pic_texture(2) = Screen('MakeTexture', handle, pic2);
n_tex = 5;
for i = 1:n_tex

eval(['load gaussians', num2str(i), '.mat']); % there are 10 different gauss files, so we load them one after the other
for ii = 1:transition_frames
    tex = round(linspace(1, 180, transition_frames));
    tex = tex(ii);
    curr_gauss = imresize(gaussians(tex).red, stimsize/500);
    pic_textures(1, ii, i) = Screen('MakeTexture', handle, cat(3, pic1, curr_gauss));
    pic_textures(2, ii, i) = Screen('MakeTexture', handle, cat(3, pic2, curr_gauss));
end
end
clear gaussians curr_gauss tex

%% Break in between runs

for i = 1:28
    message = sprintf('%i s pause', (30 - i));
    Screen('DrawText', handle, message, screen.xcen-50, screen.ycen);
    Screen('Flip', handle);
    [~, ~, keyCode, ~] = KbCheck;
    if keyCode(27)
        error('Script terminated by user');
    end
    WaitSecs(1);
end

Screen('Flip', handle);

beep;
WaitSecs(1);
beep;
WaitSecs(1);
beep;

end



%% Pracice trial: Slow transitions
for trialcount = 1:2

%% Prepare Schedule
load('event_durations.mat');
total_time = 0;
if trialcount == 1 || trialcount == 2
% values for controls
    mu_switches = 17.416;
    sd_switches = 5.01;
    mu_reversions = 1.736;
    sd_reversions = 1.19;
    DomDistr = adjDomC;
    MixDistr = adjMixC;
    meanMix = meanMixC;
elseif trialcount == 3 || trialcount == 4
    mu_switches = 17.416;
    sd_switches = 5.01;
    mu_reversions = 1.736;
    sd_reversions = 1.19;
    DomDistr = adjDomA;
    MixDistr = adjMixA;
    meanMix = meanMixA;
end


while total_time < 55 || total_time > 60
    clear transitions;
    no_of_switches = round(normrnd(mu_switches, sd_switches));
    while no_of_switches > mu_switches+sd_switches || no_of_switches < mu_switches-sd_switches
        no_of_switches = round(normrnd(mu_switches, sd_switches));
    end
    no_of_reversions = round(normrnd(mu_reversions, sd_reversions));
    while no_of_reversions > mu_reversions+sd_reversions || no_of_reversions < mu_reversions-sd_reversions
        no_of_reversions = round(normrnd(mu_reversions, sd_reversions));
    end
    
    no_of_doms = no_of_switches + no_of_reversions +1;

    dom_durations = DomDistr(randi(size(DomDistr, 1), 1, no_of_doms));
    mix_durations = MixDistr(randi(size(MixDistr, 1), 1, no_of_reversions+no_of_switches));
    transitions.type(1, 1:no_of_switches) = 1;   transitions.type(1, no_of_switches+1:no_of_switches+no_of_reversions) = 2;
    for i = 1:no_of_switches+no_of_reversions
        transitions.onset(1, i) = sum(mix_durations(1:i-1))+sum(dom_durations(1:i)); % this is the onset of this transition
        transitions.end(1, i) = sum(mix_durations(1:i))+sum(dom_durations(1:i)); % this is the end of this transition
    end

    warmup = rand*(60-transitions.end(1, no_of_switches+no_of_reversions));
    transitions.onset = transitions.onset + warmup;
    transitions.end = transitions.end + warmup;
    total_time = transitions.end(1, no_of_switches+no_of_reversions);
end
    mix_frames = round(mix_durations*60);
    transitions.type = transitions.type(randperm(no_of_switches+no_of_reversions));
    
    if trialcount == 1
    % just for the practice:
    clear transitions;
    transitions.type = [1 1 1 2 1 2];
    transitions.onset = [5 12 19 25 28 31];
    transitions.end = [8 16 23 26 30 35];
    warmup = 1;
    end
    
%% Fixed propotion (50%) at trials
j = 1;
Screen('FrameOval', handle, [0 0 0], circle_rects, 6);
Screen('Flip', handle);
KbWait;
t0 = GetSecs;

if round(rand)
    % first image will be red
    tex = randi(n_tex);
    Screen('DrawTexture', handle, pic_texture(2), [], target_rect_right);
    Screen('DrawTexture', handle, pic_texture(2), [], target_rect_left);
    Screen('DrawTexture', handle, pic_textures(1, transition_frames/2, tex), [], target_rect_right);
    Screen('DrawTexture', handle, pic_textures(1, transition_frames/2, tex), [], target_rect_left);
    Screen('FrameOval', handle, [0 0 0], circle_rects, 6);
    Screen('Flip', handle);
    
    while GetSecs < t0+warmup
    j=j+1;
    [~, secs(j, 1), keyCode, ~] = KbCheck;
    %        sacc_frames{j, 1} = GetSaccadometerNewFrame(pstr);
    keypresses(j, :) = keyCode(arrowkeys); % 38-40: arrows, 68: d, 87: w
    end
    
    Screen('DrawTexture', handle, pic_texture(1), [], target_rect_right);
    Screen('DrawTexture', handle, pic_texture(1), [], target_rect_left);
    Screen('FrameOval', handle, [0 0 0], circle_rects, 6);
    Screen('Flip', handle);
    last_image = 1;
    next_image = 2;
    
else
    % first image will be green
        % first image will be red
    tex = randi(n_tex);
    duration = rand*(60-total_time);
    Screen('DrawTexture', handle, pic_texture(1), [], target_rect_right);
    Screen('DrawTexture', handle, pic_texture(1), [], target_rect_left);
    Screen('DrawTexture', handle, pic_textures(2, transition_frames/2, tex), [], target_rect_right);
    Screen('DrawTexture', handle, pic_textures(2, transition_frames/2, tex), [], target_rect_left);
    Screen('FrameOval', handle, [0 0 0], circle_rects, 6);
    Screen('Flip', handle);
    
    while GetSecs < t0+warmup
    j=j+1;
    [~, secs(j, 1), keyCode, ~] = KbCheck;
    %        sacc_frames{j, 1} = GetSaccadometerNewFrame(pstr);
    keypresses(j, :) = keyCode(arrowkeys); % 38-40: arrows, 68: d, 87: w
    end
    
    Screen('DrawTexture', handle, pic_texture(2), [], target_rect_right);
    Screen('DrawTexture', handle, pic_texture(2), [], target_rect_left);
    Screen('FrameOval', handle, [0 0 0], circle_rects, 6);
    Screen('Flip', handle);
    last_image = 2;
    next_image = 1;

end % we start on an ambiguous image


for ii = 1:size(transitions.type,2)
    
    
    while GetSecs < t0+transitions.onset(ii)
    j=j+1;
    [~, secs(j, 1), keyCode, ~] = KbCheck;
    %        sacc_frames{j, 1} = GetSaccadometerNewFrame(pstr);
    keypresses(j, :) = keyCode(arrowkeys); % 38-40: arrows, 68: d, 87: w
    end % we wait out the dominance period

    tex = randi(n_tex);
    
    if transitions.type(ii) == 1 % switch
        
        Screen('DrawTexture', handle, pic_texture(last_image), [], target_rect_right);
        Screen('DrawTexture', handle, pic_texture(last_image), [], target_rect_left);
        Screen('DrawTexture', handle, pic_textures(next_image, transition_frames/2, tex), [], target_rect_right);
        Screen('DrawTexture', handle, pic_textures(next_image, transition_frames/2, tex), [], target_rect_left);
        Screen('FrameOval', handle, [0 0 0], circle_rects, 6);
        Screen('Flip', handle);

        while GetSecs < t0+transitions.end(ii)
        j=j+1;
        [~, secs(j, 1), keyCode, ~] = KbCheck;
        %        sacc_frames{j, 1} = GetSaccadometerNewFrame(pstr);
        keypresses(j, :) = keyCode(arrowkeys); % 38-40: arrows, 68: d, 87: w
        end % we wait out the mix period

        Screen('DrawTexture', handle, pic_texture(next_image), [], target_rect_right);
        Screen('DrawTexture', handle, pic_texture(next_image), [], target_rect_left);
        Screen('FrameOval', handle, [0 0 0], circle_rects, 6);
        Screen('Flip', handle);
        
        [next_image, last_image] = swap(next_image, last_image);
        
    elseif transitions.type(ii) == 2 % reversion
        
        Screen('DrawTexture', handle, pic_texture(last_image), [], target_rect_right);
        Screen('DrawTexture', handle, pic_texture(last_image), [], target_rect_left);
        Screen('DrawTexture', handle, pic_textures(next_image, transition_frames/2, tex), [], target_rect_right);
        Screen('DrawTexture', handle, pic_textures(next_image, transition_frames/2, tex), [], target_rect_left);
        Screen('FrameOval', handle, [0 0 0], circle_rects, 6);
        Screen('Flip', handle);

        while GetSecs < t0+transitions.end(ii)
        j=j+1;
        [~, secs(j, 1), keyCode, ~] = KbCheck;
        %        sacc_frames{j, 1} = GetSaccadometerNewFrame(pstr);
        keypresses(j, :) = keyCode(arrowkeys); % 38-40: arrows, 68: d, 87: w
        end % we wait out the mix period

        Screen('DrawTexture', handle, pic_texture(last_image), [], target_rect_right);
        Screen('DrawTexture', handle, pic_texture(last_image), [], target_rect_left);
        Screen('FrameOval', handle, [0 0 0], circle_rects, 6);
        Screen('Flip', handle);
    end
    

    
end
    
    
while GetSecs < t0+36
j=j+1;
[~, secs(j, 1), keyCode, ~] = KbCheck;
%        sacc_frames{j, 1} = GetSaccadometerNewFrame(pstr);
keypresses(j, :) = keyCode(arrowkeys); % 38-40: arrows, 68: d, 87: w
end % we wait until 60s is reached
Screen('Flip', handle);
%% Reset variables
secs = zeros(100000, 1); sacc_frames = cell(20000, 1); keypresses = zeros(100000, 3);
j = 1;

%% make new textures
pic1 = ['D:\Dropbox\Documents\University\Matlab Code\Rivalry\Simulation\Pictures\green' num2str(trialcount+1) '.JPG'];
pic2 = ['D:\Dropbox\Documents\University\Matlab Code\Rivalry\Simulation\Pictures\red' num2str(trialcount+1) '.JPG'];

pic1 = imread(pic1); pic1 = imresize(pic1, stimsize/size(pic1, 1));
pic2 = imread(pic2); pic2 = imresize(pic2, stimsize/size(pic2, 1));

pic_texture(1) = Screen('MakeTexture', handle, pic1);
pic_texture(2) = Screen('MakeTexture', handle, pic2);

n_tex = 5;
for i = 1:n_tex

eval(['load gaussians', num2str(i), '.mat']); % there are 10 different gauss files, so we load them one after the other
for ii = 1:transition_frames
    tex = round(linspace(1, 180, transition_frames));
    tex = tex(ii);
    curr_gauss = imresize(gaussians(tex).red, stimsize/500);
    pic_textures(1, ii, i) = Screen('MakeTexture', handle, cat(3, pic1, curr_gauss));
    pic_textures(2, ii, i) = Screen('MakeTexture', handle, cat(3, pic2, curr_gauss));
end
end
clear gaussians curr_gauss tex

end



%% Trial loop: RT trials
for trialcount = 1:4
j = 1;

%% Prepare Schedule
load('event_durations.mat');
total_time = 0;
if trialcount == 1 || trialcount == 2
% values for controls
    mu_switches = 17.416;
    sd_switches = 5.01;
    mu_reversions = 1.736;
    sd_reversions = 1.19;
    DomDistr = adjDomC;
    MixDistr = adjMixC;
    meanMix = meanMixC;
elseif trialcount == 3 || trialcount == 4
    mu_switches = 17.416;
    sd_switches = 5.01;
    mu_reversions = 1.736;
    sd_reversions = 1.19;
    DomDistr = adjDomA;
    MixDistr = adjMixA;
    meanMix = meanMixA;
end


while total_time < 55 || total_time > 60
    clear transitions;
    no_of_switches = round(normrnd(mu_switches, sd_switches));
    while no_of_switches > mu_switches+sd_switches || no_of_switches < mu_switches-sd_switches
        no_of_switches = round(normrnd(mu_switches, sd_switches));
    end
    no_of_reversions = round(normrnd(mu_reversions, sd_reversions));
    while no_of_reversions > mu_reversions+sd_reversions || no_of_reversions < mu_reversions-sd_reversions
        no_of_reversions = round(normrnd(mu_reversions, sd_reversions));
    end
    
    no_of_doms = no_of_switches + no_of_reversions +1;

    dom_durations = DomDistr(randi(size(DomDistr, 1), 1, no_of_doms));
    mix_durations = MixDistr(randi(size(MixDistr, 1), 1, no_of_reversions+no_of_switches));
    transitions.type(1, 1:no_of_switches) = 1;   transitions.type(1, no_of_switches+1:no_of_switches+no_of_reversions) = 2;
    for i = 1:no_of_switches+no_of_reversions
        transitions.onset(1, i) = sum(mix_durations(1:i-1))+sum(dom_durations(1:i)); % this is the onset of this transition
        transitions.end(1, i) = sum(mix_durations(1:i))+sum(dom_durations(1:i)); % this is the end of this transition
    end

    warmup = rand*(60-transitions.end(1, no_of_switches+no_of_reversions));
    transitions.onset = transitions.onset + warmup;
    transitions.end = transitions.end + warmup;
    total_time = transitions.end(1, no_of_switches+no_of_reversions);
end
mix_frames = round(mix_durations*60);
transitions.type = transitions.type(randperm(no_of_switches+no_of_reversions));

%% Fixed propotion (50%) at trials
Screen('FrameOval', handle, [0 0 0], circle_rects, 6);
Screen('Flip', handle);
KbWait;
t0 = GetSecs;

if round(rand)
    % first image will be red
    tex = randi(n_tex);
    Screen('DrawTexture', handle, pic_texture(2), [], target_rect_right);
    Screen('DrawTexture', handle, pic_texture(2), [], target_rect_left);
    Screen('DrawTexture', handle, pic_textures(1, transition_frames/2, tex), [], target_rect_right);
    Screen('DrawTexture', handle, pic_textures(1, transition_frames/2, tex), [], target_rect_left);
    Screen('FrameOval', handle, [0 0 0], circle_rects, 6);
    Screen('Flip', handle);
    
    while GetSecs < t0+warmup
    j=j+1;
    [~, secs(j, 1), keyCode, ~] = KbCheck;
    %        sacc_frames{j, 1} = GetSaccadometerNewFrame(pstr);
    keypresses(j, :) = keyCode(arrowkeys); % 38-40: arrows, 68: d, 87: w
    end
    
    Screen('DrawTexture', handle, pic_texture(1), [], target_rect_right);
    Screen('DrawTexture', handle, pic_texture(1), [], target_rect_left);
    Screen('FrameOval', handle, [0 0 0], circle_rects, 6);
    Screen('Flip', handle);
    last_image = 1;
    next_image = 2;
    first_dom = 'red';
    
else
    % first image will be green
        % first image will be red
    tex = randi(n_tex);
    duration = rand*(60-total_time);
    Screen('DrawTexture', handle, pic_texture(1), [], target_rect_right);
    Screen('DrawTexture', handle, pic_texture(1), [], target_rect_left);
    Screen('DrawTexture', handle, pic_textures(2, transition_frames/2, tex), [], target_rect_right);
    Screen('DrawTexture', handle, pic_textures(2, transition_frames/2, tex), [], target_rect_left);
    Screen('FrameOval', handle, [0 0 0], circle_rects, 6);
    Screen('Flip', handle);
    
    while GetSecs < t0+warmup
    j=j+1;
    [~, secs(j, 1), keyCode, ~] = KbCheck;
    %        sacc_frames{j, 1} = GetSaccadometerNewFrame(pstr);
    keypresses(j, :) = keyCode(arrowkeys); % 38-40: arrows, 68: d, 87: w
    end
    
    Screen('DrawTexture', handle, pic_texture(2), [], target_rect_right);
    Screen('DrawTexture', handle, pic_texture(2), [], target_rect_left);
    Screen('FrameOval', handle, [0 0 0], circle_rects, 6);
    Screen('Flip', handle);
    last_image = 2;
    next_image = 1;
    first_dom = 'green';

end % we start on an ambiguous image


for ii = 1:no_of_switches+no_of_reversions
    
    
    while GetSecs < t0+transitions.onset(ii)
    j=j+1;
    [~, secs(j, 1), keyCode, ~] = KbCheck;
    %        sacc_frames{j, 1} = GetSaccadometerNewFrame(pstr);
    keypresses(j, :) = keyCode(arrowkeys); % 38-40: arrows, 68: d, 87: w
    end % we wait out the dominance period

    tex = randi(n_tex);
    
    if transitions.type(ii) == 1 % switch
        
        Screen('DrawTexture', handle, pic_texture(last_image), [], target_rect_right);
        Screen('DrawTexture', handle, pic_texture(last_image), [], target_rect_left);
        Screen('DrawTexture', handle, pic_textures(next_image, transition_frames/2, tex), [], target_rect_right);
        Screen('DrawTexture', handle, pic_textures(next_image, transition_frames/2, tex), [], target_rect_left);
        Screen('FrameOval', handle, [0 0 0], circle_rects, 6);
        Screen('Flip', handle);

        while GetSecs < t0+transitions.end(ii)
        j=j+1;
        [~, secs(j, 1), keyCode, ~] = KbCheck;
        %        sacc_frames{j, 1} = GetSaccadometerNewFrame(pstr);
        keypresses(j, :) = keyCode(arrowkeys); % 38-40: arrows, 68: d, 87: w
        end % we wait out the mix period

        Screen('DrawTexture', handle, pic_texture(next_image), [], target_rect_right);
        Screen('DrawTexture', handle, pic_texture(next_image), [], target_rect_left);
        Screen('FrameOval', handle, [0 0 0], circle_rects, 6);
        Screen('Flip', handle);
        
        [next_image, last_image] = swap(next_image, last_image);
        
    elseif transitions.type(ii) == 2 % reversion
        
        Screen('DrawTexture', handle, pic_texture(last_image), [], target_rect_right);
        Screen('DrawTexture', handle, pic_texture(last_image), [], target_rect_left);
        Screen('DrawTexture', handle, pic_textures(next_image, transition_frames/2, tex), [], target_rect_right);
        Screen('DrawTexture', handle, pic_textures(next_image, transition_frames/2, tex), [], target_rect_left);
        Screen('FrameOval', handle, [0 0 0], circle_rects, 6);
        Screen('Flip', handle);

        while GetSecs < t0+transitions.end(ii)
        j=j+1;
        [~, secs(j, 1), keyCode, ~] = KbCheck;
        %        sacc_frames{j, 1} = GetSaccadometerNewFrame(pstr);
        keypresses(j, :) = keyCode(arrowkeys); % 38-40: arrows, 68: d, 87: w
        end % we wait out the mix period

        Screen('DrawTexture', handle, pic_texture(last_image), [], target_rect_right);
        Screen('DrawTexture', handle, pic_texture(last_image), [], target_rect_left);
        Screen('FrameOval', handle, [0 0 0], circle_rects, 6);
        Screen('Flip', handle);
    end
    

    
end
    
while GetSecs < t0+60
j=j+1;
[~, secs(j, 1), keyCode, ~] = KbCheck;
%        sacc_frames{j, 1} = GetSaccadometerNewFrame(pstr);
keypresses(j, :) = keyCode(arrowkeys); % 38-40: arrows, 68: d, 87: w
end % we wait until 60s is reached
Screen('Flip', handle);

%% Save data and reset variables
filename = strcat('D:\Dropbox\Documents\University\Matlab Code\Rivalry\Simulation\Results\', num2str(date), '_', participant_id, '_simulation_RT_trace_', num2str(trialcount), '.mat');
if iscell(filename)
    filename = filename{1};
end
secs = secs(1:j, 1); keypresses = keypresses(1:j, :); % sacc_frames = sacc_frames(1:j, 1);
secs = secs - t0;
save(filename, 'secs', 'keypresses', 'sacc_frames', 'transitions', 'mix_durations', 'dom_durations', 'first_dom');
secs = zeros(100000, 1); sacc_frames = cell(20000, 1); keypresses = zeros(100000, 3);
j = 1;

%% make new textures
pic1 = ['D:\Dropbox\Documents\University\Matlab Code\Rivalry\Simulation\Pictures\green' num2str(mod(trialcount+2, 6)+1) '.JPG'];
pic2 = ['D:\Dropbox\Documents\University\Matlab Code\Rivalry\Simulation\Pictures\red' num2str(mod(trialcount+2, 6)+1) '.JPG'];

pic1 = imread(pic1); pic1 = imresize(pic1, stimsize/size(pic1, 1));
pic2 = imread(pic2); pic2 = imresize(pic2, stimsize/size(pic2, 1));

pic_texture(1) = Screen('MakeTexture', handle, pic1);
pic_texture(2) = Screen('MakeTexture', handle, pic2);

n_tex = 5;
for i = 1:n_tex

eval(['load gaussians', num2str(i), '.mat']); % there are 10 different gauss files, so we load them one after the other
for ii = 1:transition_frames
    tex = round(linspace(1, 180, transition_frames));
    tex = tex(ii);
    curr_gauss = imresize(gaussians(tex).red, stimsize/500);
    pic_textures(1, ii, i) = Screen('MakeTexture', handle, cat(3, pic1, curr_gauss));
    pic_textures(2, ii, i) = Screen('MakeTexture', handle, cat(3, pic2, curr_gauss));
end
end
clear gaussians curr_gauss tex

%% Break in between runs
if  trialcount == 4;
    break
end
for i = 1:28
    message = sprintf('%i s pause', (30 - i));
    Screen('DrawText', handle, message, screen.xcen-d-stimsize/2, screen.ycen);
    Screen('DrawText', handle, message, screen.xcen+d-stimsize/2, screen.ycen);
    Screen('Flip', handle);
    [~, ~, keyCode, ~] = KbCheck;
    if keyCode(27)
        error('Script terminated by user');
    end
    WaitSecs(1);
end

Screen('Flip', handle);

beep;
WaitSecs(1);
beep;
WaitSecs(1);
beep;

end





%% Shutdown
Screen('CloseAll');
ShowCursor;

catch err
    ShowCursor;
    Screen('CloseAll');
    rethrow(err)
    
end