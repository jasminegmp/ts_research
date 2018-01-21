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

%% Find matrix profile
MAGIC_mp_seg_len = 100;
[matrixProfile, profileIndex, motifIdxs, discordIdx] = interactiveMatrixProfileVer2(y, MAGIC_mp_seg_len);

% Sliding window to find local minimums
local_min_value = [];
local_min_time = [];
count = 1;
for idx = 1:MAGIC_mp_seg_len:length(matrixProfile) - MAGIC_mp_seg_len
    [min_val, min_idx] = min(matrixProfile(idx:idx+MAGIC_mp_seg_len - 1));
    local_min_value(count) = min_val;
    local_min_time(count) = min_idx + idx;
    count = count + 1;
end

% Sort local minimums and its corresponding index
loc_min = transpose(vertcat(local_min_time, local_min_value));

diff_loc_min = diff(loc_min(:,1));

% for idx = 1:length(diff_loc_min)
%     if (diff_loc_min(idx) > MAGIC_mp_seg_len)
%         loc_min(idx) = [];
%     end
% end

count = 1;
exit = 1;
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
% For now, this is good enough to know where to look for distance profile

%% Debugging plots
figure;
plot(matrixProfile);
hold on;
scatter(mp_motifs(:,1),mp_motifs(:,2));
plot(y);

%% Calculate DP based on the most minimum value and consecutiveness
% That means you have to calculate the pairs of mins
[min_motif_a,min_motif_idx_a] = min(mp_motifs(:,2));

if mod(min_motif_idx_a,2) % ODD
    min_motif_idx_b = min_motif_idx_a + 1;
else % EVEN
    min_motif_idx_b = min_motif_idx_a - 1;
end
min_motif_b = mp_motifs(min_motif_idx_b,2);

% Calculate DP based on the most min value
dp = findNN(transpose(y),transpose(y(mp_motifs(min_motif_idx_a,1):mp_motifs(min_motif_idx_a,1)+MAGIC_mp_seg_len-1)));
%dp = findNN(transpose(y),transpose(y(mp_motifs(min_motif_idx_a,1):mp_motifs(min_motif_idx_a,1)+MAGIC_mp_seg_len-1)));

figure; plot(dp)
count = 1;
for idx = 1:MAGIC_mp_seg_len/2:length(dp) - MAGIC_mp_seg_len/2
    [dp_min_val, dp_min_idx] = min(dp(idx:idx+MAGIC_mp_seg_len/2 - 1));
    dp_local_min_value(count) = dp_min_val;
    dp_local_min_time(count) = dp_min_idx + idx;
    count = count + 1;
end
hold on;
scatter(dp_local_min_time, dp_local_min_value);

dp_loc_min = transpose(vertcat(dp_local_min_time, dp_local_min_value));

% Find if motif found in mp (min_motif_idx_b) exists in minimums of DP
MAGIC_threshold = .35;
threshold = round(MAGIC_threshold *MAGIC_mp_seg_len);
found = [];

for idx = 1:length(dp_loc_min)
    idx_test = abs(dp_loc_min(idx,1)-mp_motifs(min_motif_idx_b,1));
    if (idx_test <= threshold)
        found = [dp_loc_min(idx,1), dp_loc_min(idx,2)];
    end
end

%% Now you know where the motifs exist (motif_a and motif_b)
merge_a = mp_motifs(min_motif_idx_a,1);
merge_b = mp_motifs(min_motif_idx_b,1);

merged_ts = [];

merged_ts = calc_avg_ts(y, MAGIC_mp_seg_len, merge_a, merge_b);
tree = struct();

% Fill one node
[tree, merge_idx] = add_tree_node(merge_a, MAGIC_mp_seg_len,merge_idx, merged_ts, y, tree);
[tree, merge_idx] = add_tree_node(merge_b, MAGIC_mp_seg_len,merge_idx, merged_ts, y, tree);


%% Look for all other possible merges which is at dp_local_min_time and dp_local_min_value
% first remove anything in dp_local_min that is within the two merged
% values
dp_loc_min = dp_loc_min((dp_loc_min(:,1)<merge_a-MAGIC_mp_seg_len) | (dp_loc_min(:,1)>(merge_a+MAGIC_mp_seg_len)),:);
dp_loc_min = dp_loc_min((dp_loc_min(:,1)<merge_b-MAGIC_mp_seg_len) | (dp_loc_min(:,1)>(merge_b+MAGIC_mp_seg_len)),:);

% Find all minimums within MAGIC_threshold
dp_loc_min = dp_loc_min((dp_loc_min(:,2) < MAGIC_threshold*min_motif_a+min_motif_a),:);

%% Now get rid of any that are too close to each other
dp_loc_min = sort(dp_loc_min);
% Find two consecutive values, merge, delete, repeat
flag = 1;
idx = 1;
while flag
    % Find 2 consec values to merge and merge
    if (((dp_loc_min(idx+1) - dp_loc_min(idx)) > (MAGIC_mp_seg_len - MAGIC_mp_seg_len*MAGIC_threshold)) && (dp_loc_min(idx+1) - dp_loc_min(idx)) < (MAGIC_mp_seg_len + MAGIC_mp_seg_len*MAGIC_threshold))
        merged_ts = calc_avg_ts(y, MAGIC_mp_seg_len, dp_loc_min(idx,1), dp_loc_min(idx+1,1));
        [tree, merge_idx] = add_tree_node(dp_loc_min(idx,1), MAGIC_mp_seg_len,merge_idx, merged_ts, y, tree);
        [tree, merge_idx] = add_tree_node(dp_loc_min(idx+1,1), MAGIC_mp_seg_len,merge_idx, merged_ts, y, tree);
        dp_loc_min(idx,:) = [];
        dp_loc_min(idx, :) = [];
    else
        dp_loc_min(idx+1, :) = [];
    end
    if length(dp_loc_min) <= 2
        flag = 0;
    end
end

debug_plot_subtrees(tree,y)

%replace all y with parent values now


function debug_plot_subtrees(tree,y)
    figure;
    hold on;
    color = ['b','r','g','k', 'c', 'm', 'y'];
    color_idx = 1;
    for idx = 1:length(tree)
        plot((tree(idx).index(1,1):tree(idx).index(1,2)), y(tree(idx).index(1,1):tree(idx).index(1,2)), color(color_idx))
        if ~mod(idx,2) && idx > 1
            color_idx = color_idx + 1;
        end
    end
end

function merged_ts = calc_avg_ts(y, MAGIC_mp_seg_len, merge_a, merge_b)
    % Caclulate averaged data
    for idx = 0:MAGIC_mp_seg_len-1
        merged_ts(idx+1) = (y(merge_a+idx)+y(merge_b+idx))/2;
    end
end

function [tree, merge_idx] = add_tree_node(merge, MAGIC_mp_seg_len,merge_idx, merged_ts,y, tree)
    % Fill one node
    tree(merge_idx).data = y(merge:merge+MAGIC_mp_seg_len-1);
    
    tree(merge_idx).parent_idx = merged_ts;
    tree(merge_idx).index = [merge,merge+MAGIC_mp_seg_len-1];
    if merge_idx == 1 || 2
        tree(merge_idx).left_c = 0;
        tree(merge_idx).left_c = 0;
    end
    merge_idx = merge_idx + 1;
end

%%%%%%% Do this after finish merging everything
% replace time series with merged
%y(merge_a:merge_a+MAGIC_mp_seg_len) = NaN
%y(merge_b:merge_b+MAGIC_mp_seg_len) = NaN
%y(isnan(y)) = []