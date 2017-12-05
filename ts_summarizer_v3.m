%% Jasmine Kim 
% Time series summarizer - Bottom Up
% 11/13/17

%% Clean
clean_up_workspace();

%% Constant defines
segments = 8;
period = 2*pi; % Assume for now, I've already derived period from FFT
Fs = 1/period;                    % Sampling frequency
L = period*segments;                   % Length of signal
noise = 15;

%% Generate ts
[t,y] = ts_generator(Fs, L, noise);

%% Plot original time series
figure;
hold on;
plot(y);
title('Original Time Series')

%% Magic, we know subsequence length
% We know num of segments

% Break up data
sub_len = round(length(y) / segments);

% final array
summary = {}

%% create array of sub_len
ts_array = {};
% allowed for some minor overlap, hardcoded some values
for i = 1:segments
    if i == 1
        ts_array{i} = y(1: i*sub_len + 1);
    elseif i == segments
        ts_array{i} = y((i-1)*sub_len: end);
        pad_zeros = zeros(1,((i)*sub_len- (i-1)*sub_len) - (length(y)-sub_len*(i-1)));
        ts_array{i} = horzcat(ts_array{i},pad_zeros);
    else
        ts_array{i} = y((i-1)*sub_len: i*sub_len);
    end 
end
num_data_points = length(ts_array{1});

count = 1;

%% Now that we have 
level_ds = {1,segments};
merged_level_ds = {1,segments};
%segment = struct('merged_ts',[],'rearranged_ts',[])
%level_ds{1}.merged_ts = cell2struct(ts_array);
%= ts_array;
% initialize
level_ds{1} = ts_array;
merged_level_ds{1} = ts_array;
%level_ds{1,1} = ts_array{1};
merge_highlighter = [];
merge_highlighter_idx = 1

local_segments = segments;
for segment = 2:segments
    %% Create distance array
    dist_array = zeros(1,local_segments-1);
    %% Find minimum in distance array
    
    %Jasmine note - need to get this to only recalculate for updated values
    for idx = 1:length(dist_array)
        dist_array(1, idx) = findNN(ts_array{idx},ts_array{idx+1});
    end
    [min_val, min_idx] = min(dist_array);
    %% Find average of the two most similar subsequences
    for j = 1:length(ts_array{1,min_idx})
        temp{j} = (ts_array{1,min_idx}(j) + ts_array{1,min_idx+1}(j))/2;
    end
    % plot(cell2mat(temp(1,:)))
    %% Keep meta data
    original_data = transpose({ts_array{min_idx}, ts_array{min_idx+1}});
    merge_highlighter(merge_highlighter_idx) = min_idx;
    merge_highlighter_idx = merge_highlighter_idx + 1;
    %% Replace with new merged subsequence in ts_array
    ts_array{min_idx + 1} = [];
    ts_array{min_idx} = cell2mat(temp);
    ts_array = ts_array(~cellfun('isempty',ts_array));
    merged_level_ds{segment} = ts_array;

    k = segment;
    % first copy over previous level_ds into next level_ds
    level_ds{1,k} = level_ds{1,k-1};
    %remove minimum index from new cell
    temp_var = level_ds{1,k-1}{1,min_idx};
    level_ds{1,k}{1,min_idx} = [];
    level_ds{1,k} = level_ds{1,k}(~cellfun('isempty',level_ds{1,k}));
    temp_size = size(level_ds{1,k-1}{1,min_idx+1});
    temp_size_2 = size(temp_var);
    if (temp_size(1) > 1 || temp_size_2(1) > 1)
        level_ds{1,k}{1,min_idx} = vertcat(temp_var, level_ds{1,k-1}{1,min_idx+1});
    else
        level_ds{1,k}{1,min_idx} = vertcat(temp_var, {level_ds{1,k-1}{1,min_idx+1}});
    end
    
    
%end
    
    
    
   %level_ds{segment} = temp_ts_array;
    local_segments = local_segments - 1;
end

%% Now plot all of the different levels in level_ds
start = 1;
figure;
hold on;
plotNumber = 1;
merged_array = []
merge_idx = 1;
for i = 1:length(level_ds)
    subplot(segments/2, 2, plotNumber);
    axis([0 350 -2 2])
    hold on;
    start = 1;
    title(['Level ' num2str(i)]);
    for j = 1:length(level_ds{1,i})
        size_of_cell =  size(level_ds{1,i}{1,j});
        if size_of_cell(1) > 1
            for k = 1:length(level_ds{1,i}{1,j}) 
                merged_array(merge_idx) = j;
                merge_idx = merge_idx + 1;
                if ismember(k, merged_array)
                    color = [.5 .5 .5];
                    merged_array(merge_idx) = j;
                    merge_idx = merge_idx + 1;
                else
                    color = [1 0 0];
                    merged_array(merge_idx) = k;
                    merge_idx = merge_idx + 1;
                end
                if k > 1
                    %merged = [k, ]
                    plot([start: start+num_data_points-1],level_ds{1,i}{1,j}{k,1}, 'Color', color, 'LineWidth', 1);
                else
                    plot([start: start+num_data_points-1],level_ds{1,i}{1,j}{k,1}, 'Color', color, 'LineWidth', 1);
                end
            end
            start = start + num_data_points;
        else
            plot([start: start+num_data_points-1],level_ds{1,i}{1,j}, 'Color', [.5 .5 .5], 'LineWidth', 1);
            start = start + num_data_points;
        end
    end
    %merged = []
    merge_idx = 1;
    plotNumber = plotNumber + 1;
end

start = 1;
figure;
plotNumber = 1;

hold on;
found_merge = 0
merge_highlighter_idx = 1;
for j = 1:length(merged_level_ds)
    subplot(segments/2, 2, plotNumber);
    axis([0 350 -2 2])
    hold on;
    title(['Merged Level ' num2str(j)]);
    for k = 1:length(merged_level_ds{1,j})
        color = [.5 .5 .5];
        if j > 1 
            if merge_highlighter(merge_highlighter_idx) == k && found_merge == 0
                color = [1 0 0];
                found_merge = 1;
                merge_highlighter_idx = merge_highlighter_idx + 1;
            end
        end
        plot([start: start+num_data_points-1],merged_level_ds{1,j}{1,k}, 'Color', color);
        start = start + num_data_points;
    end
    found_merge = 0
    start = 1;
    plotNumber = plotNumber + 1;
end


