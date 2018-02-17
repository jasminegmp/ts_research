%function jasmine_testFastMPdist()
DEBUG = 1;

%% Clean
clean_up_workspace();
close all

%ts_1 = load('test.txt');
load('ecg_clean.mat');
ts_1 = ts;
ts_1 = transpose(ts_1);

total_length = 60000;
ts_1 = ts_1(1:total_length);
seg_len = 1000;
figure;
plot(ts_1);
%title('Time series 1');
hold on;
MAGIC_mp_seg_len = seg_len;

ts_2 = ts_1(1:seg_len-1);
plot(ts_2)
title('Time series 1 and 2');
xlim([0 length(ts_1)]);
SL = MAGIC_mp_seg_len / 2;

%% Variables before loop
merged_ts = [];

MAGIC_threshold = .40;
not_done = 1;
tree = struct();
level = 1;
distance_matrix = [];
distance_matrix_loc = [];
merge_idx = 1;
tree = struct();
merge_count = 0;
%%
% This segments the time series and calculates the distance
for idx = 1:MAGIC_mp_seg_len:total_length - MAGIC_mp_seg_len*2
    beg = idx;
    mid = idx + MAGIC_mp_seg_len;
    fin = idx + MAGIC_mp_seg_len*2;
    dist = fastMPdist_SS(ts_1(beg:mid-1), ts_1(mid:fin-1),SL);
    distance_matrix = [distance_matrix;dist];
    distance_matrix_loc = [distance_matrix_loc;beg];

end

%% Now go find minimum in distance matrix
% Merge minimum - select one of the time serieså
% update distance matrix with newly merged values - 
    % this means... update row with merged with a NaN
    % update the distance of the row BEFORE NaN and AFTER NaN
% repeat UNTIL all values in distance matrix is NaN
not_exit = 1;
merged_sum = [];

% keep a separate distance area
summation_distances = distance_matrix;

while (not_exit)
   [min_val, min_idx] = min(distance_matrix);
   % make sure next value in distance_matrix((min_idx + 1)) is not NAN
   % if it IS NAN, find keep incrementing min_idx+1 until you find one that
   % is not NAN
   next_idx = min_idx + 1;
   if isnan(distance_matrix(next_idx))
    for idx = next_idx:length(distance_matrix)
        if isnan(distance_matrix(idx))
        else
            next_idx = idx;
            break
        end
    end
    % If you get here, there's a chance entire distance matrix NaN
    if isnan(distance_matrix)
        not_exit = 0;
        break;
    end
    % If you get here, you chose last index as min_idx
    % Set min_val to NaN and continue on
    distance_matrix(min_idx) = NaN;
   end

   % Merge min_idx and next_idx
   % Add children nodes
    [tree, merge_idx] = add_tree_node(distance_matrix_loc(min_idx), MAGIC_mp_seg_len,merge_idx, ts_1, tree, level, merge_count);
    [tree, merge_idx] = add_tree_node(distance_matrix_loc(next_idx), MAGIC_mp_seg_len,merge_idx, ts_1, tree, level, merge_count);
    %add both children in seen node
    
    
   % Add merged node into tree
    merged_ts = ts_1(distance_matrix_loc(min_idx):distance_matrix_loc(min_idx)+MAGIC_mp_seg_len-1);
    [tree, merge_idx] = add_parent_node(merged_ts, distance_matrix_loc(min_idx), MAGIC_mp_seg_len,merge_idx, tree, level,merge_count);
    
    % Plot the ones I merged and highlight them.
    figure;
    title(level);
    hold on;
    plot(ts_1);
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
    
    % Now need to update distance matrix AND time series
    
    % set next_idx to Nan
    distance_matrix(next_idx) = NaN;
    
   post_idx = next_idx + 1 ; 
   if post_idx < length(distance_matrix)
       if isnan(distance_matrix(post_idx))
            for jdx = post_idx:length(distance_matrix)
                if isnan(distance_matrix(jdx))
                else
                    post_idx = jdx;
                    break
                end
            end
            if post_idx >= length(distance_matrix)
                % If you get here, there's a chance entire distance matrix NaN
                if isnan(distance_matrix)
                    not_exit = 0;
                    break;
                end
                % If you get here, you probably exhausted every posssible location
                % Set min_val to NaN and continue on
                not_exit = 0;
                break;
            end
       end

   
   % calculate new distance
     seg_1 = ts_1(distance_matrix_loc(min_idx):distance_matrix_loc(min_idx)+MAGIC_mp_seg_len-1);
     seg_1(isinf(seg_1)) = 0;
     seg_2 = ts_1(distance_matrix_loc(post_idx):distance_matrix_loc(post_idx)+MAGIC_mp_seg_len-1);
     seg_2(isinf(seg_2)) = 0;
    dist = fastMPdist_SS(seg_1, seg_2, SL);
    
    % update distance profile
    distance_matrix(min_idx) = dist(1);
    summation_distances(min_idx) = dist(1);
   end
    
   
    % remove merged from TS
    ts_1(distance_matrix_loc(next_idx):distance_matrix_loc(next_idx)+MAGIC_mp_seg_len-1) = NaN;
    ts_1(isnan(ts_1)) = Inf;    
    
    merged_sum = [merged_sum;sum(summation_distances)];
    
    % increment level
    level = level + 1;
    
end
figure;
plot(merged_sum)
title('Sum of all distances vs merge');
    

%% Functions here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function merged_ts = merge(query_1, query_2)
    merged_ts = query_1;
end
function debug_plot_subtrees(tree,y, level)
    figure;
    hold on;
    %color = ['b','r','g','k', 'c', 'm', 'y'];
    color_idx = 1;
    for idx = 1:length(tree)
        if (tree(idx).level == level) && (tree(idx).parent_idx ~= 0)
            plot((tree(idx).index(1,1):tree(idx).index(1,2)), y(tree(idx).index(1,1):tree(idx).index(1,2)));
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
%%
function [tree, merge_idx] = add_tree_node(merge, MAGIC_mp_seg_len,merge_idx,y, tree, level, merge_count)
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
    tree(merge_idx).merge_count = merge_count
end
%%
function [tree, merge_idx] = add_parent_node(merge_ts, merge, MAGIC_mp_seg_len,merge_idx, tree, level, merge_count)
    tree(merge_idx).data = merge_ts;
    tree(merge_idx).level = level + 1;
    tree(merge_idx).parent_idx = 0;
    tree(merge_idx).index = [merge,merge+MAGIC_mp_seg_len-1];
    tree(merge_idx).left_c = merge_idx-1;
    tree(merge_idx).right_c = merge_idx-2;
    tree(merge_idx).merge_count = merge_count
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

function [mp_motif_loc, mp_motif_val] = find_mp_motif_loc(loc_min, MAGIC_mp_seg_len, pattern_1_start)
    mp_motif_loc = NaN;
    mp_motif_val = NaN;
    below_range = pattern_1_start + MAGIC_mp_seg_len  - MAGIC_mp_seg_len*0.3;
    % Not checking above range because it's better to go further out
    % So that i don't run into a pathological distance of 0
   % above_range = pattern_1_start + MAGIC_mp_seg_len + MAGIC_mp_seg_len*0.8;
    for i = 1:length(loc_min)
        if (loc_min(i,1) > below_range)
            mp_motif_loc = loc_min(i,1);
            mp_motif_val = loc_min(i,2);
            break;
        end
    end
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