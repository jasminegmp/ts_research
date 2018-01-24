%% Jasmine Kim 
% Time series summarizer - Bottom up using matrix profile
% 1/10/2018

%% Clean
clean_up_workspace();

%% Constant defines
segments = 8;
period = 2*pi; % Assume for now, I've already derived period from FFT
Fs = 1/period;                    % Sampling frequency
L = period*segments*10;                   % Length of signal
noise = 4;
merge_idx = 1;
 
%% Generate ts
[t,y] = ts_generator(Fs, L, noise);

%% Plot original time series
figure;
hold on;
plot(y);
title('Original Time Series')

%% Variables before loop
merged_ts = [];
MAGIC_mp_seg_len = 100;
MAGIC_threshold = .40;
not_done = 1;
tree = struct();
level = 1;
unmodified_y = y;
%% Loop begins HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while not_done
    %% Find matrix profile
    [matrixProfile, profileIndex, motifIdxs, discordIdx] = interactiveMatrixProfileVer2(y, MAGIC_mp_seg_len);

    %% Find minimums in matrix profile

    % Sliding window to find local minimums
    loc_min = find_local_minimums(matrixProfile, MAGIC_mp_seg_len);

    % Find mp motif locations to test against DP
    mp_motifs = find_mp_motif_loc(loc_min, MAGIC_mp_seg_len);

    %% Debugging plots
    figure;
    plot(matrixProfile);
    hold on;
    scatter(mp_motifs(:,1),mp_motifs(:,2));
    plot(y);


    %% Now we can compare found motif and DP
    % Calculate DP based on the most minimum value and consecutiveness
    % That means you have to calculate the pairs of mins
    
    % If motif size is only 1, you know you can't reduce anymore. Done!
    if length(mp_motifs) == 2
        break;
    end
    
    [min_motif_a,min_motif_idx_a, min_motif_b, min_motif_idx_b] = find_min_motif(mp_motifs);

    % Calculate DP based on the most min value
    dp = findNN(transpose(y),transpose(y(mp_motifs(min_motif_idx_a,1):mp_motifs(min_motif_idx_a,1)+MAGIC_mp_seg_len-1)));

    % Plot DP
    figure; plot(dp);

    % Find all minimums in DP
    dp_loc_min = find_dp_min(MAGIC_mp_seg_len, dp);

    % Find if motif found in mp (min_motif_idx_b) exists in minimums of DP
    found = motif_in_dp(MAGIC_threshold, MAGIC_mp_seg_len, mp_motifs, min_motif_idx_b, dp_loc_min);
%%
    if ~isempty(found)
        %% Now you know where first merge exists (motif_a and motif_b)
        merge_a = mp_motifs(min_motif_idx_a,1);
        merge_b = mp_motifs(min_motif_idx_b,1);

        merged_ts = calc_avg_ts(y, MAGIC_mp_seg_len, merge_a, merge_b);
        merged_start = min(merge_a, merge_b);


        % Fill one node
        [tree, merge_idx] = add_tree_node(merge_a, MAGIC_mp_seg_len,merge_idx, y, tree, level);
        [tree, merge_idx] = add_tree_node(merge_b, MAGIC_mp_seg_len,merge_idx, y, tree, level);
        [tree, merge_idx] = add_parent_node(merged_ts, merged_start, MAGIC_mp_seg_len,merge_idx, tree, level);

        %% Look for all other possible merges from DP
        % first remove anything in dp_local_min that is within the two merged
        % values
        dp_loc_min = dp_loc_min((dp_loc_min(:,1)<merge_a-MAGIC_mp_seg_len) | (dp_loc_min(:,1)>(merge_a+MAGIC_mp_seg_len)),:);
        dp_loc_min = dp_loc_min((dp_loc_min(:,1)<merge_b-MAGIC_mp_seg_len) | (dp_loc_min(:,1)>(merge_b+MAGIC_mp_seg_len)),:);

        % Find all minimums within MAGIC_threshold
        dp_loc_min = dp_loc_min((dp_loc_min(:,2) <= MAGIC_threshold*min_motif_a+min_motif_a),:);



        %% Now recursively merge the ones that are within DP right now

        % Now get rid of any minimums in DP that are too close to each other AND merge

        %dp_loc_min = clean_dp(dp_loc_min, MAGIC_mp_seg_len, MAGIC_threshold);


        flag = 1;
        while flag
            % make sure dp_loc_min isn't empty, otherwise exist
            if isempty(dp_loc_min)
                break;
            end
            if length(dp_loc_min) == 2 % only has one value
                break;
            end
            %look for values to merge
            found = 0;
            for idx = 1:length(dp_loc_min)- 1
                if dp_loc_min(idx + 1, 1) - dp_loc_min(idx,1) >= MAGIC_mp_seg_len - MAGIC_mp_seg_len*MAGIC_threshold
                    if dp_loc_min(idx + 1, 1) - dp_loc_min(idx,1) <= MAGIC_mp_seg_len + MAGIC_mp_seg_len*MAGIC_threshold
                        found = idx;
                        break;
                    end
                end
            end
            if found == 0 % didn't find any value to merge, so exit
                break;
            end
            % found value to merge
            merged_ts = calc_avg_ts(y, MAGIC_mp_seg_len, dp_loc_min(found,1), dp_loc_min(found+1,1));
            merged_start = min(dp_loc_min(found,1), dp_loc_min(found+1,1));
            [tree, merge_idx] = add_tree_node(dp_loc_min(found,1), MAGIC_mp_seg_len,merge_idx, y, tree, level);
            [tree, merge_idx] = add_tree_node(dp_loc_min(found+1,1), MAGIC_mp_seg_len,merge_idx, y, tree, level);
            [tree, merge_idx] = add_parent_node(merged_ts, merged_start, MAGIC_mp_seg_len,merge_idx, tree, level);
            % pop off merged values
            dp_loc_min(found,:) = [];
            dp_loc_min(found, :) = [];
            % repeat
        end



        %% Replace y with merged
        figure;
        hold on;
        plot(y, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.2);
        for idx = 1:length(tree)
            if (tree(idx).parent_idx == 0) && (tree(idx).level == level + 1)
                y(tree(tree(idx).left_c).index(1,1):tree(tree(idx).left_c).index(1,2)) = NaN;
                y(tree(tree(idx).right_c).index(1,1):tree(tree(idx).left_c).index(1,2)) = NaN;
                replace_idx = min(tree(tree(idx).right_c).index(1,1), tree(tree(idx).right_c).index(1,1));
                y(replace_idx:replace_idx+length(tree(idx).data)-1) = tree(idx).data;
                y(isnan(y)) = [];
            end
        end
        
        %plot(unmodified_y);
        title(level);
        for idx = 1:length(tree)
            if (tree(idx).level == level) || (tree(idx).level == level+1)
                time = tree(idx).index(1,1):1:tree(idx).index(1,2);
                if (tree(idx).parent_idx == 0) && (tree(idx).level == level+1)
                    plot(time, tree(idx).data, 'Color', [1 0 0], 'LineWidth', 1);
                elseif (tree(idx).level == level) && (tree(idx).parent_idx ~= 0)
                    plot(time, tree(idx).data, 'Color', [0 0 1], 'LineWidth', 0.7);
                end
            end
        end
    end
    level = level + 1;
end

% LOOP END


%% Functions here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function debug_plot_subtrees(tree,y, level)
    figure;
    hold on;
    %color = ['b','r','g','k', 'c', 'm', 'y'];
    color_idx = 1;
    for idx = 1:length(tree)
        if (tree(idx).level == level) && (tree(idx).parent_idx ~= 0)
            plot((tree(idx).index(1,1):tree(idx).index(1,2)), y(tree(idx).index(1,1):tree(idx).index(1,2)))
            if ~mod(idx,2) && idx > 1
                color_idx = color_idx + 1;
            end
        end
    end
end

function merged_ts = calc_avg_ts(y, MAGIC_mp_seg_len, merge_a, merge_b)
    % Caclulate averaged data
    for idx = 0:MAGIC_mp_seg_len-1
        merged_ts(idx+1) = (y(merge_a+idx)+y(merge_b+idx))/2;
    end
end

function [tree, merge_idx] = add_tree_node(merge, MAGIC_mp_seg_len,merge_idx,y, tree, level)
    % Fill one node
    tree(merge_idx).data = y(merge:merge+MAGIC_mp_seg_len-1);
    tree(merge_idx).level = level;
    if (mod(merge_idx,3) == 1) 
        tree(merge_idx).parent_idx = merge_idx + 2;
    else % EVEN
        tree(merge_idx).parent_idx = merge_idx + 1;
    end

    tree(merge_idx).index = [merge,merge+MAGIC_mp_seg_len-1];
    if merge_idx == 1 || 2
        tree(merge_idx).left_c = 0;
        tree(merge_idx).right_c = 0;
    end
    merge_idx = merge_idx + 1;
end

function [tree, merge_idx] = add_parent_node(merge_ts, merge, MAGIC_mp_seg_len,merge_idx, tree, level)
    tree(merge_idx).data = merge_ts;
    tree(merge_idx).level = level + 1;
    tree(merge_idx).parent_idx = 0;
    tree(merge_idx).index = [merge,merge+MAGIC_mp_seg_len-1];
    tree(merge_idx).left_c = merge_idx-1;
    tree(merge_idx).right_c = merge_idx-2;
    merge_idx = merge_idx + 1;
end

function loc_min = find_local_minimums(ts, MAGIC_mp_seg_len)
    % Sliding window to find local minimums
    local_min_value = [];
    local_min_time = [];
    count_idx = 1;
    for loop_idx = 1:MAGIC_mp_seg_len:length(ts) - MAGIC_mp_seg_len
        [loc_min_val, loc_min_idx] = min(ts(loop_idx:loop_idx+MAGIC_mp_seg_len - 1));
        local_min_value(count_idx) = loc_min_val;
        local_min_time(count_idx) = loc_min_idx + loop_idx;
        count_idx = count_idx + 1;
    end
    % Sort local minimums and its corresponding index
    loc_min = transpose(vertcat(local_min_time, local_min_value));
end

function mp_motifs = find_mp_motif_loc(loc_min, MAGIC_mp_seg_len)
    count = 1;
    exit = 1;
    diff_loc_min = diff(loc_min(:,1));
    while(exit)
        if ((diff_loc_min(count) > MAGIC_mp_seg_len) || (diff_loc_min(count) < MAGIC_mp_seg_len/2))
            loc_min(count,:) = [];
            diff_loc_min = diff(loc_min(:,1));
        else
            count = count + 2;
        end
        if count >= length(diff_loc_min)
            exit = 0;
        end
    end

    mp_motifs = sortrows(loc_min, 1, 'ascend');
end

function [min_motif_a,min_motif_idx_a, min_motif_b, min_motif_idx_b] = find_min_motif(mp_motifs)
    [min_motif_a,min_motif_idx_a] = min(mp_motifs(:,2));
    if mod(min_motif_idx_a,2) % ODD
        min_motif_idx_b = min_motif_idx_a + 1;
    else % EVEN
        min_motif_idx_b = min_motif_idx_a - 1;
    end
    min_motif_b = mp_motifs(min_motif_idx_b,2);
end

function dp_loc_min = find_dp_min(MAGIC_mp_seg_len, dp)
    count = 1;
    for idx = 1:MAGIC_mp_seg_len/2:length(dp) - MAGIC_mp_seg_len/2
        [dp_min_val, dp_min_idx] = min(dp(idx:idx+MAGIC_mp_seg_len/2 - 1));
        dp_local_min_value(count) = dp_min_val;
        dp_local_min_time(count) = dp_min_idx + idx;
        count = count + 1;
    end
    dp_loc_min = transpose(vertcat(dp_local_min_time, dp_local_min_value));
end

function found = motif_in_dp(MAGIC_threshold, MAGIC_mp_seg_len, mp_motifs, min_motif_idx_b, dp_loc_min)
    % Find if motif found in mp (min_motif_idx_b) exists in minimums of DP
    threshold = round(MAGIC_threshold *MAGIC_mp_seg_len);
    found = [];
    for idx = 1:length(dp_loc_min)
        idx_test = abs(dp_loc_min(idx,1)-mp_motifs(min_motif_idx_b,1));
        if (idx_test <= threshold)
            found = [dp_loc_min(idx,1), dp_loc_min(idx,2)];
        end
    end
end

function dp_loc_min = clean_dp(dp_loc_min, MAGIC_mp_seg_len, MAGIC_threshold)
    diff_dp = diff(dp_loc_min);
    for idx = 1:2:length(diff_dp)/2
        if (diff_dp(idx,1) >= MAGIC_mp_seg_len - MAGIC_mp_seg_len*MAGIC_threshold) && (diff_dp(idx,1) <= MAGIC_mp_seg_len + MAGIC_mp_seg_len*MAGIC_threshold)
        else
            dp_loc_min(idx,:) = NaN;
            dp_loc_min(idx+1, :) = NaN;
        end
    end
    dp_loc_min(isnan(dp_loc_min)) = [];
    
end



%%%%%%% Do this after finish merging everything
% replace time series with merged
%y(merge_a:merge_a+MAGIC_mp_seg_len) = NaN
%y(merge_b:merge_b+MAGIC_mp_seg_len) = NaN
%y(isnan(y)) = []