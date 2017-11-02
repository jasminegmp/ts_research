%% Jasmine Kim 
% Time series summarizer - A.1-4 working
% 10/26/17

%% Clean
clean_up_workspace();

%% Constant defines
period = 2*pi; % Assume for now, I've already derived period from FFT
tot_subplot = 8;
curr_subplot = 1;
Fs = 10;                    % Sampling frequency
L = 1000;                   % Length of signal
noise = 1;

%% Generate ts
[t,y] = ts_generator(Fs, L, noise);

%% Prep for MASS
y = transpose(y);
m = Fs*round(period); % Subsequence length

%% Plot original time series
figure;
subplot(tot_subplot,1,curr_subplot);
curr_subplot = curr_subplot + 1;
hold on;
plot(y);
title('Original Time Series')

%discretieze ts


%% Calculate and plot Matrix Profile of time series
figure(1);
subplot(tot_subplot,1,curr_subplot);
MP = Time_series_Self_Join_Fast(y,m);
curr_subplot = curr_subplot + 1;
hold on;
plot(MP);
title('Matrix Profile')

%% Begin Algorithm
flag = 1;
curr_subseq = 1;
end_index = 1;
while(flag)
    %%
    

    %% Calculate distance profile of subsequence
     if end_index > 1
         y = y(end_index:end);
     end
    DP = findNN(y,y(1:m));
    figure(1)
    subplot(tot_subplot,1,curr_subplot);
    hold on;
    plot(DP);
    H_plot = gca;
    H_plot.XLim = [0 1000];
    title('Distance Profile')
    
    %% Find snap-on x-axis locations
    TF = islocalmin(DP,'MinProminence',2);
    % When calculating Minima, only keep the minima that is within 20% of
    % each other
    Minima = DP(TF);
    MinIdx = find(TF);
    remove_indices = find(Minima(1)*(1.2) < Minima);
    Minima(remove_indices) = [];
    MinIdx(remove_indices) = [];
    length(MinIdx)
    plot(MinIdx, Minima,'o');
    curr_subplot = curr_subplot + 1;
    
    %% Distance measure between pattern and query
    figure;
    matching_count = 1;
    pattern = y(1:MinIdx(1));
    hold on;
    plot(pattern, 'Color','r','LineWidth',2);
    H_plot_2 = gca;
    H_plot_2.XLim = [0 100];
    H_plot_2.YLim = [-5 5];
    dist_query_pattern = [];
    query_index = [];
    query_cell = {};
    for subseq_count = 1:length(MinIdx)
        query_temp = y(MinIdx(subseq_count):MinIdx(subseq_count) + length(pattern) - 1);
        query_cell{subseq_count} = query_temp;
        dist_query_pattern(subseq_count) = findNN(pattern, query_cell{subseq_count});
    end
    match_flag = 1;
    for dist_query_count = 1:length(dist_query_pattern)
        if (dist_query_pattern(dist_query_count) < dist_query_pattern(1)*(1.2))
            plot_handle = plot(query_cell{dist_query_count}, 'color', [.5 .5 .5], 'LineWidth', 1);
            plot_handle.Color(4) = 0.5;
            matching_count = matching_count + 1;
            end_index = MinIdx(dist_query_count) + m;
        end
    end
    plot(pattern, 'Color','r','LineWidth',2);
    
    title(['Similar Finds: ' num2str(matching_count)]);
    
    %% Finished going through time series
    if length(y) < end_index + m
        flag = 0;
    end
    
end
