%% Clean
clear
clc
delete(findall(0,'Type','figure'))
close all

%% Modifiable constants
DEBUG = 1;
default_file = '12726m.mat';
segment_length = 1000;

%% Input file
% expecting all load files to have a ts time series datapoint

if DEBUG
    data = load(default_file);
    ts = data.val;
    ts = transpose(ts);
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

fastMPdist_seg_len = round(seg_len / 2);

dist_mat = {};
removed_seg = [];
merge_info = struct();
merged_tree = struct();
ts_1_nan = ts_1;

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
count = 1;
%% Divide up time series by segment length and claculate distance array
for idx = 1:seg_len:tot_len - seg_len*2
    
    seg_0 = idx;
    seg_1 = idx + seg_len - 1;
    seg_2 = seg_1 + seg_len - 1;

    dist = fastMPdist_SS(ts_1(seg_0:seg_1), ts_1(seg_1:seg_2),fastMPdist_seg_len);
    seg_ts_1 = [ts_1(seg_0:seg_1)];
    seg_ts_2 = [ts_1(seg_1:seg_2)];
    dist_info = [dist seg_0 seg_1 seg_2];
    dist_mat{count,1} = dist;
    dist_mat{count,2} = seg_0;
    dist_mat{count,3} = seg_1;
    dist_mat{count,4} = seg_2;
    dist_mat{count,5} = seg_ts_1;
    dist_mat{count,6} = seg_ts_2;
    count = count + 1;
    %dist_mat = {dist_mat;dist_info};
    
end

%% keep track of merge meta data
merge_info.start = cell2mat(dist_mat(:,2));
merge_info.end = cell2mat(dist_mat(:,3));
merge_info.merge_count = zeros(length(dist_mat),1);


%temp var
final_dist = dist_mat;
merge_count = 1;
seen_dist = [];
dist_mat(:,7) = dist_mat(:,1);
%% Now this is the main loop that....
% 1. Goes find minimum in distance matrix
% 2. Merges minimum
% 3. Updates distance matrix
% 4. Repeat until all time series segments have been merged
while (size(dist_mat,1) > 1)
    
    % 1.
    [min_val, min_idx] = min(cell2mat(dist_mat(:,1)));
    min_idx = min_idx(1);
    min_val = min_val(1);
    

    % 2.
    % found segments to be merged
    loc_0 = cell2mat(dist_mat(min_idx,2));
    loc_1 = cell2mat(dist_mat(min_idx,3));
    loc_2 = cell2mat(dist_mat(min_idx,4));
    m_seg_0 = ts_1(loc_0:loc_1);
    m_seg_1 = ts_1(loc_1:loc_2);
%     
%     if isnan(current) && isnan(previous)
%         current = loc_1;
%         previous = current; 
%         plot_loc_0 = loc_0; 
%         plot_loc_1 = loc_1;
%         plot_loc_2 = loc_2;
%     else
%         if previous >= current
%             plot_loc_0 = loc_0; 
%            % plot_loc_0 = loc_0 + segment_length*realign_count; 
%             plot_loc_1 = loc_1 + segment_length*realign_count;
%             plot_loc_2 = loc_2 + segment_length*realign_count;
%             realign_count = realign_count + 1;
%         end
%         current = loc_1;
%         previous = current; 
%     end
    
    
%     % plot the two segments that are going to be merged
%     figure;
%     title(level);
%     hold on;
%     plot(ts_1,'Color', [.1 .4 1], 'LineWidth', 0.7);
%     plot(loc_0:loc_1, m_seg_0, 'Color', [1 0 0], 'LineWidth', 0.7)
%     plot(loc_1:loc_2, m_seg_1, 'Color', [0 0 1], 'LineWidth', 0.7)


    % merge one segment into another
    %ts_1(loc_1:loc_2) = [];
    ts_1(loc_1:loc_2) = [];
    
    % meta data info
    
    index  = strfind(transpose(ts), transpose(m_seg_1));
    merged_idx = find(merge_info.start == (index+1));
    merge_info.merge_count(merged_idx) = level;
    
    % place new info into tree struct
    merged_tree = add_tree_node(merged_tree, level, min_val, cell2mat(dist_mat(min_idx,2)), cell2mat(dist_mat(min_idx,3)), cell2mat(dist_mat(min_idx,4)), m_seg_0, m_seg_1, ts_1, sum(cell2mat(dist_mat(:,1))));
    level = level + 1;

    % 3.
    % replace distance at min_dix to NaN so it's not used again
    dist_mat(min_idx,:) = [];
    diff_mat = setdiff(cell2mat(final_dist(:,1)), cell2mat(dist_mat(:,7)));
    diff_val = setdiff(diff_mat, seen_dist);
    seen_dist = vertcat(diff_val,seen_dist);
    seen_dist = unique(seen_dist);
    row_mark = find(cell2mat(final_dist(:,1)) == diff_val);
    final_dist(row_mark,7) = {merge_count};
    merge_count = merge_count + 1;
    if min_idx <= size(dist_mat,1)
        % update time values
        for idx = min_idx:size(dist_mat,1)
            for j = 2:4 %updating columns 2 -4
                dist_mat(idx,j) = {cell2mat(dist_mat(idx,j)) - seg_len};
            end
        end
        % update distance at min_idx + 1 unless it's the last item in dist_mat
        dist = fastMPdist_SS(ts_1(cell2mat(dist_mat(min_idx,2)):cell2mat(dist_mat(min_idx,3))), ts_1(cell2mat(dist_mat(min_idx,3)):cell2mat(dist_mat(min_idx,4))),fastMPdist_seg_len);
        dist_mat(min_idx,1) = {dist};
    end
    % 4.
    
end

% % plot final two segments that are going to be merged
% loc_0 = dist_mat(1,2);
% loc_1 = dist_mat(1,3);
% loc_2 = dist_mat(1,4);
% m_seg_0 = ts_1(loc_0:loc_1);
% m_seg_1 = ts_1(loc_1:loc_2);

% figure;
% title(level);
% hold on;
% plot(ts_1);
% plot(loc_0:loc_1, m_seg_0, 'Color', [0 0 1], 'LineWidth', 0.7)
% plot(loc_1:loc_2, m_seg_1, 'Color', [1 0 0], 'LineWidth', 0.7)
% ts_1 = ts_1(loc_0:loc_1);
% merged_tree = add_tree_node(merged_tree, level, min_val, loc_0, loc_1, loc_2, m_seg_0, m_seg_1, ts_1, dist_mat(1,1));
% level = level + 1;

% Go plot output
for j = 1:length(final_dist)
    figure;
    title(j);
    hold on;
    plot(ts,'Color', [.1 .4 1], 'LineWidth', 0.7);
    
    for k = 1:1:length(final_dist)
        if cell2mat(final_dist(k,7)) <= j
            plot(cell2mat(final_dist(k,2)):cell2mat(final_dist(k,3)), cell2mat(final_dist(k,5)), 'Color', [0.5 0.5 0.5], 'LineWidth', 0.7)
        end
    end
    
    %plot(loc_1:loc_2, m_seg_1, 'Color', [0 0 1], 'LineWidth', 0.7)
end
    


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