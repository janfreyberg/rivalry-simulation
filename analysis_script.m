t=1;
s=1;
files = dir('*.mat');

prev_subj = '';
icRR = 1;
icRG = 1;
icGG = 1;
icGR = 1;
irRR = 1;
irRG = 1;
irGG = 1;
irGR = 1;

%%
for i = 1:size(files, 1)
clear switch_time switch_type

        load(files(i).name);

    curr_file = strsplit(files(i).name, '_');
    subject = curr_file{2};
    trial_type = curr_file{4};
    trial_no = str2num(curr_file{6}(1));

    if ~strcmp(subject, prev_subj)
        icRR = 1;
        icRG = 1;
        icGG = 1;
        icGR = 1;
        irRR = 1;
        irRG = 1;
        irGG = 1;
        irGR = 1;
        prev_subj = subject;
    end


%%
    observations = size(keypresses, 1);
    switches = size(transitions.type, 2);
    j=1;


    for ii = 2:observations
        if ~isequal(keypresses(ii, :), keypresses(ii-1, :)) && ~isequal(keypresses(ii, :), [0 0 0])
            switch_type(j, :) = keypresses(ii, :);
            switch_time(j) = secs(ii);
            j=j+1;
        end
    end

    if strcmp(first_dom, 'green')
    last_dominance = 3;
    elseif strcmp(first_dom, 'red')
    last_dominance = 1;
    end


%%
    for ii = 1:switches-1
        
        if last_dominance == 1 && transitions.type(ii) == 1
            next_dominance = 3;
        elseif last_dominance == 3 && transitions.type(ii) == 1 
            next_dominance = 1;
        elseif last_dominance == 1 && transitions.type(ii) == 2
            next_dominance = 1;
        elseif last_dominance == 3 && transitions.type(ii) == 2
            next_dominance = 3;
        end
        
        
        onset_index = find(switch_time > transitions.onset(ii), 1);
        ons_betw = 0;
        while ~switch_type(onset_index, 2) || switch_time(onset_index)-transitions.onset(ii) < 0.250
            
            if ~switch_type(onset_index, last_dominance) || (sum(switch_type(onset_index, :))==2 && ~switch_type(onset_index, 2))
            ons_betw = ons_betw + 1;
            disp('      ii=   last_dominance=');
            disp([ii last_dominance]);
            disp(switch_type(onset_index, :));
            end
            onset_index = onset_index+1;
        end
        
        end_index = onset_index+1;
        end_betw = 0;
        if strcmp(trial_type, 'RT')
        end_index = find(switch_time > transitions.end(ii), 1);
        end
        
            if end_index > size(switch_type, 1)
                end_index = size(switch_type, 1);
            end
        while ~switch_type(end_index, next_dominance) || isequal(switch_type(end_index, :), [1 0 1])
            if end_index > size(switch_type, 1)
                end_index = size(switch_type, 1);
                break
            end
            
            if switch_type(end_index, :) ~= [0 1 0]
            end_betw = end_betw + 1;
            disp('      ii=   next_dominance=');
            disp([ii next_dominance]);
            disp(switch_type(end_index, :));
            end
            
            end_index = end_index+1;
            
            if end_index > size(switch_type, 1)
                end_index = size(switch_type, 1);
                break
            end
        end
        
        
        rx_1 = switch_time(onset_index) - transitions.onset(ii);
        rx_2 = switch_time(end_index) - transitions.onset(ii);
        rx_3 = switch_time(end_index) - transitions.end(ii);
        
        
        if switch_time(onset_index) > transitions.onset(ii+1);
            ons_missed = 1;
        else
            ons_missed = 0;
        end
        if switch_time(end_index) > transitions.end(ii+1);
            end_missed = 1;
        else
            end_missed = 0;
        end
        
        s=s+1;
        
        if strcmp(trial_type, 'CRITERION')
        if trial_no == 1 || trial_no == 2
            rx_1 = rx_1/1.4060;
            rx_2 = rx_2/1.4060;
        elseif trial_no == 3 || trial_no == 4
            rx_1 = rx_1/2.1517;
            rx_2 = rx_2/2.1517;
        end
        end
        
%%        
if strcmp(trial_type, 'RT')
        if last_dominance == 1 && next_dominance == 3
            RT.switch_rt_GA(irGR) = rx_1;
            RT.switch_rt_AR(irGR) = rx_3;
            RT.switch_missed_GA(irGR) = ons_missed;
            RT.switch_missed_AR(irGR) = end_missed;
            irGR = irGR+1;
        elseif last_dominance == 3 && next_dominance == 1
            RT.switch_rt_RA(irRG) = rx_1;
            RT.switch_rt_AG(irRG) = rx_3;
            RT.switch_missed_RA(irRG) = ons_missed;
            RT.switch_missed_AG(irRG) = end_missed;
            irRG = irRG+1;
        elseif last_dominance == 1 && next_dominance == 1
            RT.rev_rt_GA(irGG) = rx_1;
            RT.rev_rt_AG(irGG) = rx_3;
            RT.rev_missed_GA(irGG) = ons_missed;
            RT.rev_missed_AG(irGG) = end_missed;
            irGG = irGG+1;
        elseif last_dominance == 3 && next_dominance == 3
            RT.rev_rt_RA(irRR) = rx_1;
            RT.rev_rt_AR(irRR) = rx_3;
            RT.rev_missed_RA(irRR) = ons_missed;
            RT.rev_missed_AR(irRR) = end_missed;
            irRR = irRR+1;
        end
elseif strcmp(trial_type, 'CRITERION')
        if last_dominance == 1 && next_dominance == 3
            CRIT.switch_rt_GA(icGR) = rx_1;
            CRIT.switch_rt_AR(icGR) = rx_2;
            CRIT.switch_missed_GA(icGR) = ons_missed;
            CRIT.switch_missed_AR(icGR) = end_missed;
            icGR = icGR+1;
        elseif last_dominance == 3 && next_dominance == 1
            CRIT.switch_rt_RA(icRG) = rx_1;
            CRIT.switch_rt_AG(icRG) = rx_2;
            CRIT.switch_missed_RA(icRG) = ons_missed;
            CRIT.switch_missed_AG(icRG) = end_missed;
            icRG = icRG+1;
        elseif last_dominance == 1 && next_dominance == 1
            CRIT.rev_rt_GA(icGG) = rx_1;
            CRIT.rev_rt_AG(icGG) = rx_2;
            CRIT.rev_missed_GA(icGG) = ons_missed;
            CRIT.rev_missed_AG(icGG) = end_missed;
            icGG = icGG+1;
        elseif last_dominance == 3 && next_dominance == 3
            CRIT.rev_rt_RA(icRR) = rx_1;
            CRIT.rev_rt_AR(icRR) = rx_2;
            CRIT.rev_missed_RA(icRR) = ons_missed;
            CRIT.rev_missed_AR(icRR) = end_missed;
            icRR = icRR+1;
        end
end


if strcmp(trial_type, 'RT')
        if last_dominance == 1 && next_dominance == 3
            RT.switch_betw_GA(irGR) = ons_betw;
            RT.switch_betw_AR(irGR) = end_betw;
            irGR = irGR+1;
        elseif last_dominance == 3 && next_dominance == 1
            RT.switch_betw_RA(irRG) = ons_betw;
            RT.switch_betw_AG(irRG) = end_betw;
            irRG = irRG+1;
        elseif last_dominance == 1 && next_dominance == 1
            RT.rev_betw_GA(irGG) = ons_betw;
            RT.rev_betw_AG(irGG) = end_betw;
            irGG = irGG+1;
        elseif last_dominance == 3 && next_dominance == 3
            RT.rev_betw_RA(irRR) = ons_betw;
            RT.rev_betw_AR(irRR) = end_betw;
            irRR = irRR+1;
        end
elseif strcmp(trial_type, 'CRITERION')
        if last_dominance == 1 && next_dominance == 3
            CRIT.switch_betw_GA(icGR) = ons_betw;
            CRIT.switch_betw_AR(icGR) = end_betw;
            icGR = icGR+1;
        elseif last_dominance == 3 && next_dominance == 1
            CRIT.switch_betw_RA(icRG) = ons_betw;
            CRIT.switch_betw_AG(icRG) = end_betw;
            icRG = icRG+1;
        elseif last_dominance == 1 && next_dominance == 1
            CRIT.rev_betw_GA(icGG) = ons_betw;
            CRIT.rev_betw_AG(icGG) = end_betw;
            icGG = icGG+1;
        elseif last_dominance == 3 && next_dominance == 3
            CRIT.rev_betw_RA(icRR) = ons_betw;
            CRIT.rev_betw_AR(icRR) = end_betw;
            icRR = icRR+1;
        end
end








if last_dominance ~= next_dominance
    last_dominance = next_dominance;
end

%%






        
    end
if strcmp(trial_type, 'RT')
eval([subject '_' trial_type '= RT;']);
elseif strcmp(trial_type, 'CRITERION')
eval([subject '_' trial_type '= CRIT;']);
end

end