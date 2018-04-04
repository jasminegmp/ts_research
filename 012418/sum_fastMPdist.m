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
needs_merging = 1;

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
    seg_2 = idx + seg_len*2 - 2;

    dist = fastMPdist_SS(ts_1(seg_0:seg_1), ts_1(seg_1:seg_2),fastMPdist_seg_len);
    dist_info = [dist seg_0 seg_1 seg_2];
    dist_mat = [dist_mat;dist_info];
    
end

%% Now this is the main loop that....
% 1. Goes find minimum in distance matrix
% 2. Merges minimum
% 3. Updates distance matrix
% 4. Repeat until all time series segments have been merged
while (needs_merging)
    [min_val, min_idx] = min(dist_mat);
    
end
_