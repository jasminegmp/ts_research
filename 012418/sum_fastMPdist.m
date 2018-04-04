%% Clean
clear
clc
delete(findall(0,'Type','figure'))
close all

%% Modifiable constants
DEBUG = 1;
default_file = 'ecg_clean.mat';
segment_length = 1000;

%% Input file
% expecting all load files to have a ts time series datapoint

if DEBUG
    data = load(default_file);
    ts = data.ts;
    seg_len = segment_length;
else
    prompt = 'Input file name: ';
    input_file = input(prompt,'s');
    prompt_2 = 'Input segment length. Note! Must be divisible by the total length. ';
    seg_len = input(prompt_2);
    
    data = load(input_file);
    ts = data.ts;
end

%% Initialization
tot_len = length(ts);
ts_1 = ts;
ts_2 = ts_1(1:seg_len);

fastMPdist_seg_len = seg_len / 2;

dist_mat = [];
merged_tree = struct();

level = 0;
needs_merging = 1;
merged_tree = add_tree_node(merged_tree, level, NaN, NaN, NaN, NaN, NaN, NaN, ts_1, NaN);
level = level + 1;

%% Plot initial time series and segment
figure;
hold on;
plot(ts_1);
plot(ts_2);
title('Time series and Segment');

%% Divide up time series by segment length and claculate distance array
for idx = 1:seg_len:tot_len - seg_len*2
    seg_0 = idx;
    seg_1 = idx + seg_len - 1;
    seg_2 = seg_1 + seg_len - 1;

    dist = fastMPdist_SS(ts_1(seg_0:seg_1), ts_1(seg_1:seg_2),fastMPdist_seg_len);
    dist_info = [dist seg_0 seg_1 seg_2];
    dist_mat = [dist_mat;dist_info];
    
end

%% Now this is the main loop that....
% 1. Goes find minimum in distance matrix
% 2. Merges minimum
% 3. Updates distance matrix
% 4. Repeat until all time series segments have been merged
while (size(dist_mat,1) > 1)
    
    % 1.
    [min_val, min_idx] = min(dist_mat);
    min_idx = min_idx(1);
    min_val = min_val(1);
    
    % 2.
    % found segments to be merged
    loc_0 = dist_mat(min_idx,2);
    loc_1 = dist_mat(min_idx,3);
    loc_2 = dist_mat(min_idx,4);
    m_seg_0 = ts_1(loc_0:loc_1);
    m_seg_1 = ts_1(loc_1:loc_2);
    
    % plot the two segments that are going to be merged
    figure;
    title(level);
    hold on;
    plot(ts_1);
    plot(loc_0:loc_1, m_seg_0, 'Color', [0 0 1], 'LineWidth', 0.7)
    plot(loc_1:loc_2, m_seg_1, 'Color', [1 0 0], 'LineWidth', 0.7)
    
    % merge one segment into another
    ts_1(loc_1:loc_2) = [];
    % place new info into tree struct
    merged_tree = add_tree_node(merged_tree, level, min_val, dist_mat(min_idx,2), dist_mat(min_idx,3), dist_mat(min_idx,4), m_seg_0, m_seg_1, ts_1, sum(dist_mat(:,1)));
    level = level + 1;
    
    % 3.
    % replace distance at min_dix to NaN so it's not used again
    dist_mat(min_idx,:) = [];
    if min_idx <= size(dist_mat,1)
        % update time values
        for idx = min_idx:size(dist_mat,1)
            for j = 2:4 %updating columns 2 -4
                dist_mat(idx,j) = dist_mat(idx,j) - 1000;
            end
        end
        % update distance at min_idx + 1 unless it's the last item in dist_mat
        dist = fastMPdist_SS(ts_1(dist_mat(min_idx,2):dist_mat(min_idx,3)), ts_1(dist_mat(min_idx,3):dist_mat(min_idx,4)),fastMPdist_seg_len);
        dist_mat(min_idx,1) = dist;
    end
    % 4.
    
end

% plot final two segments that are going to be merged
loc_0 = dist_mat(1,2);
loc_1 = dist_mat(1,3);
loc_2 = dist_mat(1,4);
m_seg_0 = ts_1(loc_0:loc_1);
m_seg_1 = ts_1(loc_1:loc_2);

figure;
title(level);
hold on;
plot(ts_1);
plot(loc_0:loc_1, m_seg_0, 'Color', [0 0 1], 'LineWidth', 0.7)
plot(loc_1:loc_2, m_seg_1, 'Color', [1 0 0], 'LineWidth', 0.7)
ts_1 = ts_1(loc_0:loc_1);
merged_tree = add_tree_node(merged_tree, level, min_val, loc_0, loc_1, loc_2, m_seg_0, m_seg_1, ts_1, dist_mat(1,1));
level = level + 1;

%plot final 
figure;
title(level);
hold on;
plot(ts_1)


% Fill one node
function tree = add_tree_node(tree, level, dist, merge_pt_0,merge_pt_1,merge_pt_2, m_seg_0, m_seg_1, ts, dist_sum)

    merge_idx = level + 1;
    
    tree(merge_idx).level = level;
    tree(merge_idx).dist = dist;
    tree(merge_idx).merge_pt_0 = merge_pt_0;
    tree(merge_idx).merge_pt_1 = merge_pt_1;
    tree(merge_idx).merge_pt_2 = merge_pt_2;
    tree(merge_idx).m_seg_0 = m_seg_0;
    tree(merge_idx).m_seg_1 = m_seg_1;
    tree(merge_idx).ts = ts;
    tree(merge_idx).dist_sum = dist_sum;
    
    
end