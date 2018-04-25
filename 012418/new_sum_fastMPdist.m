%% Clean
clear
clc
delete(findall(0,'Type','figure'))
close all

%% Modifiable constants
DEBUG = 1;
%default_file = '12726m.mat';
%segment_length = 1000;
merge_type = 'left_linkage'; % 'left_linkage', 'average_linkage', 'minimum_linkage'
fid = fopen('test3.txt');
A=textscan(fid,'%s');
ts = [];
for a = 1:length(A{1,1})
    ts(a) = str2num((A{1,1}{a,1}));
end
transpose(ts);
segment_length = 199;
seg_len = segment_length;

%% Input file
% expecting all load files to have a ts time series datapoint

% if DEBUG
%     data = load(default_file);
%     ts = data.val;
%     ts = transpose(ts);
%     seg_len = segment_length;
% else
%     prompt = 'Input file name: ';
%     input_file = input(prompt,'s');
%     prompt_2 = 'Input segment length. Note! Must be divisible by the total length. ';
%     seg_len = input(prompt_2);
%     
%     data = load(input_file);
%     ts = data.ts;
% end

%% Initialization
tot_len = length(ts);
ts_1 = ts;
ts_2 = ts_1(1:seg_len);

fastMPdist_seg_len = round(seg_len / 2);

dist_mat = {};
merge_history = {};
merged_data = {};
removed_seg = [];
merge_info = struct();
merged_tree = struct();
ts_1_nan = ts_1;

level = 0;
needs_merging = 1;
%merged_tree = add_tree_node(merged_tree, level, NaN, NaN, NaN, NaN, NaN, NaN, ts_1, NaN);
level = level + 1;
merge_count = 1;

%% Plot initial time series and segment
figure;
hold on;
plot(ts_1);
plot(ts_2);
title('Time series and Segment');
count = 1;
%% Divide up time series by segment length and claculate distance array
for idx = 1:seg_len:tot_len - seg_len
    
    seg_0 = idx;
    seg_1 = idx + seg_len;
    seg_2 = seg_1;
    seg_3 = seg_1 + seg_len;

    if seg_3 > tot_len
        break;
    end
    dist = fastMPdist_SS(ts_1(seg_0:seg_1), ts_1(seg_2:seg_3),fastMPdist_seg_len);
    seg_ts_1 = [ts_1(seg_0:seg_1)];
    seg_ts_2 = [ts_1(seg_2:seg_3)];
    dist_info = [dist seg_0 seg_1 seg_2 seg_3];
    dist_mat{count,1} = dist;
    dist_mat{count,2} = seg_0;
    dist_mat{count,3} = seg_1;
    dist_mat{count,4} = seg_2;
    dist_mat{count,5} = seg_3;
    %dist_mat{count,6} = seg_ts_1;
    %dist_mat{count,7} = seg_ts_2;
    count = count + 1;
    %dist_mat = {dist_mat;dist_info};
    
end

% create giant matrix with all merge history
merge_history{merge_count,1} = dist_mat;


%% Now this is the main loop that....
% 1. Goes find minimum in distance matrix
% 2. Merges minimum
% 3. Updates distance matrix
% 4. Repeat until all time series segments have been merged
while (size(dist_mat,1) > 1)
    if any(~isnan(cell2mat(dist_mat(:,1)))) == 0
        break;
    end

    last_element = 0;
    
    % 1. find minimum in distance matrix
    [min_val, min_idx] = min(cell2mat(dist_mat(:,1)));
    min_idx = min_idx(1);
    min_val = min_val(1);
    

    % find segments to be merged
    loc_0 = cell2mat(dist_mat(min_idx,2));
    loc_1 = cell2mat(dist_mat(min_idx,3));
    loc_2 = cell2mat(dist_mat(min_idx,4));
    loc_3 = cell2mat(dist_mat(min_idx,5));
    m_seg_0 = ts_1(loc_0:loc_1);
    m_seg_1 = ts_1(loc_2:loc_3);
    merge_history{merge_count, 1} = dist_mat;
 
   

    % check what type of merge is going to be used for next set of data
    
    %%
    % LEFT LINKAGE
    if strcmp('left_linkage',merge_type)
        
        %plot the two segments that are going to be merged
        figure;
        title(merge_count);
        hold on;
        plot(ts_1, 'LineWidth', 0.7);
        plot(loc_0:loc_1, m_seg_0, 'Color', [1 0 0], 'LineWidth', 1.2)
        plot(loc_2:loc_3, m_seg_1, 'Color', [0.1 0.1 0.1], 'LineWidth', 1.2)

        %plot old as gray
        if merge_count > 1
            for i = 1:size(merged_data,1)
                temp_loc_0 = cell2mat(merged_data{i, 1});
                temp_loc_1 = cell2mat(merged_data{i, 2});
                temp_seg = cell2mat(merged_data{i, 3});
                plot(temp_loc_0:temp_loc_1, temp_seg, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.7)
            end
        end
        
        %update distance matrix with new next values
        if min_idx < length(dist_mat)
            % find next valid value
            for j = min_idx+1:length(dist_mat)
                if ~isnan(cell2mat(dist_mat(j,1)))
                    next_loc_0 = cell2mat(dist_mat(j,4));
                    next_loc_1 = cell2mat(dist_mat(j,5));
                    new_m_seg_1 = ts_1(next_loc_0:next_loc_1);
                    merge_loc = j;
                    if min_idx+1 == j
                        last_element = 0;
                    else
                        dist_mat(j,2) = dist_mat(min_idx,2);
                        dist_mat(j,3) = dist_mat(min_idx,3);
                        last_element = 0;
                    end
                    break;
                else
                    last_element = 1;
                end
            end
            % if only nan left
        end
    
        % 2.  perform merge
        % ts_1(loc_0:loc_1)
        merged_data{merge_count,1} = {loc_2};
        merged_data{merge_count,2} = {loc_3};
        merged_data{merge_count,3} = {ts_1(loc_2:loc_3)};
        ts_1(loc_2+1:loc_3-1) = NaN;
        
        % 3. update distance matrix, calculate new distance and location
        % Remove merged segment
        dist_mat(min_idx,:) = {NaN};
        dist_mat{min_idx,8} = merge_count;
        
        if min_idx < length(dist_mat) && (last_element == 0)
            % remove all NaN and replace with 0
            m_seg_0(isnan(m_seg_0)) = 0;
            m_seg_1(isnan(m_seg_1)) = 0;
            dist = fastMPdist_SS(m_seg_0,new_m_seg_1,fastMPdist_seg_len);
            dist_mat(merge_loc,1) = {dist};
            dist_mat(merge_loc,2) = {loc_0};
            dist_mat(merge_loc,3) = {loc_1};
        end 
    end
    %%
    % AVERAGE LINKAGE
    if strcmp('average_linkage',merge_type)
        
        %update distance matrix with new next values
        if min_idx < length(dist_mat)
            %find next valid value
            for j = min_idx+1:length(dist_mat)
                if ~isnan(cell2mat(dist_mat(j,1)))
                    next_loc_0 = cell2mat(dist_mat(j,4));
                    next_loc_1 = cell2mat(dist_mat(j,5));
                    new_m_seg_1 = ts_1(next_loc_0:next_loc_1);
                    merge_loc = j;
                    if min_idx+1 == j
                        last_element = 0;
                    else
                        dist_mat(j,2) = dist_mat(min_idx,2);
                        dist_mat(j,3) = dist_mat(min_idx,3);
                        last_element = 0;
                    end
                    break;
                else
                    last_element = 1;
                end
            end
            %if only nan left
        end
    
        
        new_ts = [];
        % Find new average subsequence
        for k = 1:segment_length+1
            new_ts(k) = (ts_1(loc_0+k-1) +ts_1(loc_2+k-1))/2;
        end
        
        % 2.  perform merge, placing new averaged subsequence into earlier
        % time
        merged_data{merge_count,1} = {loc_0};
        merged_data{merge_count,2} = {loc_1};
        merged_data{merge_count,3} = {ts_1(loc_0:loc_1)};
        merged_data{merge_count,4} = {loc_2};
        merged_data{merge_count,5} = {loc_3};
        merged_data{merge_count,6} = {ts_1(loc_2:loc_3)};
        merged_data{merge_count,7} = {loc_0};
        merged_data{merge_count,8} = {loc_1};
        merged_data{merge_count,9} = {new_ts};
        merged_data{merge_count,10} = {min_val};
        ts_1(loc_0:loc_1) = new_ts;
        ts_1(loc_2:loc_3) = NaN;
        
        % 3. update distance matrix, calculate new distance and location
        % Remove merged segment
        dist_mat(min_idx,:) = {NaN};
        dist_mat{min_idx,8} = merge_count;

        if min_idx < length(dist_mat) && (last_element == 0)
            % remove all NaN and replace with 0
            m_seg_0 = merged_data{merge_count,9}{1,1};
            %m_seg_1 = merged_data{merge_count,9}{1,1};
            m_seg_0(isnan(m_seg_0)) = 0;
           % m_seg_1(isnan(m_seg_1)) = 0;
            dist = fastMPdist_SS(m_seg_0,new_m_seg_1,fastMPdist_seg_len);
            dist_mat(merge_loc,1) = {dist};
            dist_mat(merge_loc,2) = {loc_0};
            dist_mat(merge_loc,3) = {loc_1};
        end
        
        % 4. plot everything that happened in this step
        
        %plot the two segments that are going to be merged
        figure;
        title(merge_count);
        hold on;
        plot(ts_1, 'LineWidth', 0.7);
        plot(merged_data{merge_count,1}{1,1}:merged_data{merge_count,2}{1,1}, merged_data{merge_count,3}{1,1}, 'Color', [0.1 0.1 0.1], 'LineWidth', 0.7)
        plot(merged_data{merge_count,4}{1,1}:merged_data{merge_count,5}{1,1}, merged_data{merge_count,6}{1,1}, 'Color', [0.1 0.1 0.1], 'LineWidth', 0.7)
        plot(merged_data{merge_count,7}{1,1}:merged_data{merge_count,8}{1,1}, merged_data{merge_count,9}{1,1}, 'Color', [1 0 0], 'LineWidth', 1)
        %plot old as gray
        if merge_count > 1
            for i = 1:size(merged_data,1)
                temp_loc_0 = cell2mat(merged_data{i, 4});
                temp_loc_1 = cell2mat(merged_data{i, 5});
                temp_seg = cell2mat(merged_data{i, 6});
                plot(temp_loc_0:temp_loc_1, temp_seg, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.7)
            end
        end

    end
    %%
    if strcmp('minimum_linkage',merge_type)
        
        % ts_1(loc_0:loc_1)
        merged_data{merge_count,1} = {loc_2};
        merged_data{merge_count,2} = {loc_3};
        merged_data{merge_count,3} = {ts_1(loc_2:loc_3)};
        ts_1(loc_2:loc_3) = NaN;
        
        % Remove merged segment
        dist_mat(min_idx,:) = {NaN};
        dist_mat{min_idx,8} = merge_count;

        if min_idx < length(dist_mat) && (last_element == 0)
        %3. calculate new distance and location
            % remove all NaN and replace with 0
            m_seg_0(isnan(m_seg_0)) = 0;
            m_seg_1(isnan(m_seg_1)) = 0;
            dist = fastMPdist_SS(m_seg_0,new_m_seg_1,fastMPdist_seg_len);
            dist_mat(merge_loc,1) = {dist};
            dist_mat(merge_loc,2) = {loc_0};
            dist_mat(merge_loc,3) = {loc_1};
        end
    
    end
    
    merge_count = merge_count + 1;

    
end
