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

figure; plot(dp)

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
MAGIC_threshold = .1
if exists()
end

% If it does, now go find all minimums that is within threshold of the
% mp_motif
    % Summarize time series now
    % For all the minimums found ex: [48, 100, 158, 204]
    % Go and merge these and build a merged array
    % Remove merged array from original time series now
% go and re run matrix profile and re-generate a matrix profile










