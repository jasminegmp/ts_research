% The prototype for faster distance calculation
% Shaghayegh Gharghabi / Eamonn Keogh 08/29/2017
%
% [distance] = ...
%     fasterMPdist_SS(Ts1, Ts2, subLen, Thr);
% Output:
%     distance: the distance between Ts1 and Ts2 with different length(scalar)
% Input:
%     Ts1: the first input time series (vector)
%     Ts2: the second input time series (vector)
%     subLen: subsequence length (scalar)
%     Thr: threshol for distance in MP
% This code does not use STAMP or STOMP, it uses SCRIMP

function distance = fastMPdist_SS( Ts1, Ts2 , SL)
% Ts1 the longer time series
% Ts2 the longer time series
 Thr = 0.05;
tic
MASS_dist = calcMassDistMatrix(Ts1, Ts2, SL);
[MASSDist_row, MASSDist_col] = size(MASS_dist);

% calculate row-min and column-min
MASS_distSlidMin(1:MASSDist_row, 1:MASSDist_col) = 0;
allRightHistogram = min(MASS_dist);
for i = 1:MASSDist_row
    MASS_distSlidMin(i,:) = movmin(MASS_dist(i,:), MASSDist_row);%length(Ts2)-SL+1 ~MASSDist_row
end

MPdistLength = length(Ts1)-length(Ts2)+1; rightHistLength = length(Ts2)-SL+1; %% length of left and right hist is equal
MPdistArray(1:MPdistLength-1) = 0;
leftHist(1:rightHistLength) = 0; %rightHist(1:rightHistLength) = 0;

for i = 1 : MPdistLength
    rightHist = allRightHistogram(i:rightHistLength+i-1); %% maybe you can write it faster
    leftHist(:) = MASS_distSlidMin(:, i+ceil(MASSDist_row/2)); %min(MASS_dist(j,(i+length(Ts2)-1)-SL+1));
    %--------------- calc distance in specific location -----------------%
    recreatedMP = [leftHist, rightHist];
    MPdistArray(i) = calMPdist_withoutInf(recreatedMP, Thr, 2*length(Ts2));
end

toc
%figure
%plot(MPdistArray);
%title('MPdist using MASS');
distance = MPdistArray;
end

function MASS_dist = calcMassDistMatrix(Ts1, Ts2, SL)

numofSubSeq = length(Ts2)-SL +1;
MASS_dist(1:numofSubSeq, 1:length(Ts1)-SL+1) = 0; %initialize
for i = 1:numofSubSeq
    MASS_dist(i,:) = real(MASS_V4(Ts1, Ts2(i:i+SL-1))); %MASS_V4(Ts1', Ts2(i:i+SL-1)')
end

end