
%% Clean
clear
clc
delete(findall(0,'Type','figure'))
close all

%% Modifiable constants
DEBUG = 1;
%default_file = '12726m.mat';
%segment_length = 1000;
merge_type = 'average_linkage'; % 'left_linkage', 'average_linkage', 'minimum_linkage'
fid = fopen('test2.txt');
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
s_plot_count = 1;
orig_ts = ts_1;
figure;
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
        if s_plot_count == 5
            figure;
            s_plot_count = 1;
        end
        
        subplot(4,1,s_plot_count);
        s_plot_count = s_plot_count + 1;
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
    
                %find prev value
        if min_idx > 1
            prev = NaN;
            %find next valid value
            for j = min_idx-1:-1:1
                if ~isnan(cell2mat(dist_mat(j,1)))
                    if min_idx+1 == j
                        last_element = 0;
                    else
                        last_element = 0;
                    end
                    prev = j;
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
        old_ts = ts_1;
        ts_1(loc_0:loc_1) = new_ts;
        ts_1(loc_2:loc_3) = NaN;
        
        % 3. update distance matrix, calculate new distance and location
        % Remove merged segment
        dist_mat(min_idx,:) = {NaN};
        dist_mat{min_idx,8} = merge_count;

        m_seg_0 = merged_data{merge_count,9}{1,1};
        %m_seg_1 = merged_data{merge_count,9}{1,1};
        m_seg_0(isnan(m_seg_0)) = 0;
        if min_idx < length(dist_mat)
            % remove all NaN and replace with 0

           % m_seg_1(isnan(m_seg_1)) = 0;
            dist = fastMPdist_SS(m_seg_0,new_m_seg_1,fastMPdist_seg_len);
            dist_mat(merge_loc,1) = {dist};
            dist_mat(merge_loc,2) = {loc_0};
            dist_mat(merge_loc,3) = {loc_1};
            
        end
        if min_idx > 1 && ~isnan(prev)
                        
            prev_seg_loc_0 = cell2mat(dist_mat(prev,2));
            prev_seg_loc_1 = cell2mat(dist_mat(prev,3));
            prev_seg_0 = ts(prev_seg_loc_0:prev_seg_loc_1);
            prev_seg_0(isnan(prev_seg_0)) = 0;
            dist = fastMPdist_SS(m_seg_0,prev_seg_0,fastMPdist_seg_len);
            dist_mat(prev,1) = {dist};
            
        end
        
        
        % 4. plot everything that happened in this step
        
        %plot the two segments that are going to be merged
        if s_plot_count == 5
            figure;
            s_plot_count = 1;
        end
        
        subplot(4,1,s_plot_count);
        s_plot_count = s_plot_count + 1;
        
        %plot the two segments that are going to be merged
        %figure;
        title(merge_count);
        hold on;
        
        plot(orig_ts, 'LineWidth', 0.5, 'LineWidth', 1.5,'Color', [.8 .8 .65]);
        plot(old_ts, 'LineWidth', 0.7, 'Color',[0.25 0.42 .73]);
        plot(merged_data{merge_count,1}{1,1}:merged_data{merge_count,2}{1,1}, merged_data{merge_count,3}{1,1}, 'Color', [0.4 0.4 0.4], 'LineWidth', 0.7)
        plot(merged_data{merge_count,4}{1,1}:merged_data{merge_count,5}{1,1}, merged_data{merge_count,6}{1,1}, 'Color', [0.4 0.4 0.4], 'LineWidth', 0.7)
        plot(merged_data{merge_count,7}{1,1}:merged_data{merge_count,8}{1,1}, merged_data{merge_count,9}{1,1}, 'Color', [1 0 0], 'LineWidth', 1)
%         %plot old as gray
%         if merge_count > 1
%             for i = 1:size(merged_data,1)
%                 if i < merge_count
%                     temp_loc_0 = cell2mat(merged_data{i, 4});
%                     temp_loc_1 = cell2mat(merged_data{i, 5});
%                     temp_seg = cell2mat(merged_data{i, 6});
%                     plot(temp_loc_0:temp_loc_1, temp_seg, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.7)
%                 end
%             end
%         end

    end
    %%
    if strcmp('minimum_linkage',merge_type)
        
        %update distance matrix with new next values
        if min_idx < length(dist_mat)
            next = NaN;
            %find next valid value
            for j = min_idx+1:length(dist_mat)
                if ~isnan(cell2mat(dist_mat(j,1)))
                    if min_idx+1 == j
                        last_element = 0;
                    else
                        last_element = 0;
                    end
                    next = j;
                    break;
                else
                    last_element = 1;
                    
                end
            end
            %if only nan left
        end
        
        %find prev value
        if min_idx > 1
            prev = NaN;
            %find next valid value
            for j = min_idx-1:-1:1
                if ~isnan(cell2mat(dist_mat(j,1)))
                    if min_idx+1 == j
                        last_element = 0;
                    else
                        last_element = 0;
                    end
                    prev = j;
                    break;
                else
                    last_element = 1;
                    
                end
            end
            %if only nan left
        end
    

        
        % 2.  perform merge, replacing the greater distance subsequence 
        % old
        merged_data{merge_count,1} = {loc_0};
        merged_data{merge_count,2} = {loc_1};
        merged_data{merge_count,3} = {ts_1(loc_0:loc_1)};
        merged_data{merge_count,4} = {loc_2};
        merged_data{merge_count,5} = {loc_3};
        merged_data{merge_count,6} = {ts_1(loc_2:loc_3)};
        
        old_ts = ts_1;
        
        if isnan(prev) && isnan(next)
            break;
        end
        
        % find which one to replace
        if min_idx <= length(dist_mat)
            % after
            if next <= length(dist_mat)
                current_m_seg_1 = cell2mat({ts_1(loc_2:loc_3)});
                temp_m_seg_1 = cell2mat({ts_1(cell2mat(dist_mat(next,4)):cell2mat(dist_mat(next,5)))});
                current_m_seg_1(isnan(current_m_seg_1)) = 0;
                temp_m_seg_1(isnan(temp_m_seg_1)) = 0;
            end
                
            if min_idx > 1 && ~isnan(prev)
                % before
                current_m_seg_0 = cell2mat({ts_1(loc_0:loc_1)});
                temp_m_seg_0 = cell2mat({ts_1(cell2mat(dist_mat(prev,2)):cell2mat(dist_mat(prev,3)))});
                current_m_seg_0(isnan(current_m_seg_0)) = 0;
                temp_m_seg_0(isnan(temp_m_seg_0)) = 0;
            
            end

            % find if we want to merge before or after
            if min_idx > 1 && next > min_idx && ~isnan(prev)
                % compare
                dist_0 = fastMPdist_SS(current_m_seg_0,temp_m_seg_0,fastMPdist_seg_len);
                dist_1 = fastMPdist_SS(current_m_seg_1,temp_m_seg_1,fastMPdist_seg_len);
                
                
                % the after is better
                if dist_0 > dist_1
                    m_idx = next;
                    r_idx = prev; %removed idx
                    loc_a = loc_2;
                    loc_b = loc_3;
                    loc_c = loc_0;
                    loc_d = loc_1;
                    dist = dist_1;
                % the before is better
                else
                    m_idx = prev;
                    r_idx = next; %removed idx
                    loc_a = loc_0;
                    loc_b = loc_1;
                    loc_c = loc_2;
                    loc_d = loc_3;
                    dist = dist_0;
                end   
            % has no after  
            elseif min_idx == next && ~isnan(prev)
                m_idx = prev;
                r_idx = min_idx;
                dist = fastMPdist_SS(current_m_seg_0,temp_m_seg_0,fastMPdist_seg_len);
                loc_a = loc_0;
                loc_b = loc_1;
                loc_c = loc_2;
                loc_d = loc_3;
                dist = dist_0;
            % HAS NO BEFORE     
            else
                m_idx = next;
                r_idx = min_idx;
                dist = fastMPdist_SS(current_m_seg_1,temp_m_seg_1,fastMPdist_seg_len);
                loc_a = loc_2;
                loc_b = loc_3;
                loc_c = loc_0;
                loc_d = loc_1;
                dist = dist_1;
                
            end


            
            % calculate new distance
            if next ~= prev && next > prev
                if m_idx == next
                    dist_mat(r_idx,4) = {loc_a};
                    dist_mat(r_idx,5) = {loc_b};
                else
                    dist_mat(r_idx,2) = {loc_a};
                    dist_mat(r_idx,3) = {loc_b};

                end
            end
            
            
            m_seg_0 = cell2mat({ts_1(cell2mat(dist_mat(m_idx,2)):cell2mat(dist_mat(m_idx,3)))});
            m_seg_1 = cell2mat({ts_1(cell2mat(dist_mat(m_idx,4)):cell2mat(dist_mat(m_idx,5)))});
            m_seg_0(isnan(m_seg_0)) = 0;
            m_seg_1(isnan(m_seg_1)) = 0;
            dist = fastMPdist_SS(m_seg_0,m_seg_1,fastMPdist_seg_len);
            
            dist_mat(m_idx,1) = {dist};
            m_ts = ts_1(loc_c:loc_d);
            new_ts = ts_1(loc_a:loc_b);
            ts_1(loc_c:loc_d) = NaN;
            
            
        end
        


        

        
        %new
        merged_data{merge_count,7} = {loc_a};
        merged_data{merge_count,8} = {loc_b};
        merged_data{merge_count,9} = {new_ts};
        merged_data{merge_count,10} = {min_val};
        merged_data{merge_count,11} = {loc_c};
        merged_data{merge_count,12} = {loc_d};
        merged_data{merge_count,13} = {m_ts};

        
        
        
        
        % 3. update distance matrix, calculate new distance and location
        % Remove merged segment
        dist_mat(min_idx,:) = {NaN};
        dist_mat{min_idx,8} = merge_count;


        
        % 4. plot everything that happened in this step
        
        %plot the two segments that are going to be merged
        if s_plot_count == 5
            figure;
            s_plot_count = 1;
        end
        
        subplot(4,1,s_plot_count);
        s_plot_count = s_plot_count + 1;
        
        %plot the two segments that are going to be merged
        %figure;
        title(merge_count);
        hold on;
        
       % plot(orig_ts, 'LineWidth', 0.5, 'LineWidth', 1.5,'Color', [.8 .8 .65]);
        plot(ts_1, 'LineWidth', 0.7, 'Color',[0.25 0.42 .73]);
        plot(merged_data{merge_count,1}{1,1}:merged_data{merge_count,2}{1,1}, merged_data{merge_count,3}{1,1}, 'Color', [0.2 0.2 0.2], 'LineWidth', 0.7)
        plot(merged_data{merge_count,4}{1,1}:merged_data{merge_count,5}{1,1}, merged_data{merge_count,6}{1,1}, 'Color', [0.2 0.2 0.2], 'LineWidth', 0.7)
        plot(merged_data{merge_count,7}{1,1}:merged_data{merge_count,8}{1,1}, merged_data{merge_count,9}{1,1}, 'Color', [1 0 0], 'LineWidth', 1)
    
        %plot old as gray
        if merge_count > 1
            for i = 1:size(merged_data,1)
                if i < merge_count
                    temp_loc_0 = cell2mat(merged_data{i, 11});
                    temp_loc_1 = cell2mat(merged_data{i, 12});
                    temp_seg = cell2mat(merged_data{i, 13});
                    plot(temp_loc_0:temp_loc_1, temp_seg, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.7)
                end
            end
        end
        
    end
    
    merge_count = merge_count + 1;

    
end
