% The prototype for distance calculation
% Shaghayegh Gharghabi / Eamonn Keogh 08/29/2017
%
% [distance] = ...
%     MP_distance(Ts1, Ts2, subLen, Thr);
% Output:
%     distance: the distance between Ts1 and Ts2 with same length(scalar)
% Input:
%     Ts1: the first input time series (vector)
%     Ts2: the second input time series (vector)
%     subLen: subsequence length (scalar)
%     Thr: threshol for distance in MP
% This code does not use STAMP or STOMP, it uses SCRIMP

function distance = MP_distance( Ts1, Ts2 , SubLen, Thr)
if (length(Ts1) ~= length(Ts2))
    error(['Error: The length of two times series are not the same']);
end

[row, ~] = size(Ts1);
if row ~=1
    Ts = [Ts1',Ts2'];
else
    Ts = [Ts1,Ts2];
end

changePoint = length(Ts1);
[MP, ~] = MatrixProfileSplitConstraint(Ts, SubLen, changePoint);

% Calc MP Dist
MPLength = length(Ts);
distLoc = ceil(Thr*MPLength);
MPSorted = sort(MP);

if MPSorted(distLoc)~= inf
    distance = MPSorted(distLoc);
else
    MPRemoveInf = MP(MP(:)~=inf);
    distance = max(MPRemoveInf);
end   
end
