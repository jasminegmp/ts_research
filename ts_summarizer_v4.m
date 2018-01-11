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
noise = 1;
 
%% Generate ts
[t,y] = ts_generator(Fs, L, noise);

%% Plot original time series
figure;
hold on;
plot(y);
title('Original Time Series')

%% Find matrix profile
MAGIC_mp_seg_len = 50;
[matrixProfile, profileIndex, motifIdxs, discordIdx] = interactiveMatrixProfileVer2(y, MAGIC_mp_seg_len);

%% Find the top 10% minimum matrix profile values
MAGIC_percentage = 0.1;
n = round(MAGIC_percentage*length(y));
[mp_y,mp_t] = mink(matrixProfile, n, 1);

%% Go through minimums looking for first dips that are within MP segment length
%  These minimum values are stored in target_dp
target_dp = [];
diff_mp_t = diff(mp_t);
for idx = 1:diff_mp_t
    if abs(diff_mp_t(idx)) <= MAGIC_mp_seg_len
        target_dp = [target_dp, idx];
    end
end

%% 
dp = findNN(y,y(target_dp(1):target_dp(1)+MAGIC_mp_seg_len));