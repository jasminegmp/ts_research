% The prototype for distance calculation
% Shaghayegh Gharghabi / Eamonn Keogh 08/29/2017
%
% [distance] = ...
%     MP_distance_SS(Ts1, Ts2, subLen, Thr);
% Output:
%     distance: the distance between Ts1 and Ts2 with different length(scalar)
% Input:
%     Ts1: the first input time series (vector)
%     Ts2: the second input time series (vector)
%     subLen: subsequence length (scalar)
%     Thr: threshol for distance in MP
% This code does not use STAMP or STOMP, it uses SCRIMP

function distance = MPdist_SS( Ts1, Ts2 , SubLen, Thr)

l1 = length(Ts1); l2 = length(Ts2);
tic
distance(1:abs(l1-l2+1)) = 0;
if l1 < l2
    for i = 1:l2-l1+1
        distance(i) = MPdist(Ts1, Ts2(i: i+l1-1), SubLen, Thr);
    end
else
    for i = 1:l1-l2+1
        distance(i) = MPdist(Ts2, Ts1(i: i+l2-1), SubLen, Thr);

    end
end
toc


% xlimit = max(l1 , l2);
figure
% subplot(3,1,1);
% plot(Ts1);
% title('TS1');
% xlim([0, xlimit]);
% subplot(3,1,2);
% plot(Ts2);
% title('TS2');
% xlim([0, xlimit]);
% subplot(3,1,3);
plot(distance);
title('MPdist using MP');
% xlim([0, xlimit]);

end