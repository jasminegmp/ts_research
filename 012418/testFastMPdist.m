function testFastMPdist()
Ts1 = load('test.txt');
close all

Ts1 = Ts1(1:1000);
Ts2 = Ts1(500:500+100); %awgn(-x', 10) -x'
SL = 10;
figure
plot(Ts1);
title('Ts1');

figure
plot(Ts2);
title('Ts2');
xlim([0 length(Ts1)]);

distanceFastMPdist = fastMPdist_SS(Ts1, Ts2, SL);
% distanceMPB = MPdist_SS(Ts1, Ts2 , SL, Thr);
% difference = diff(distanceFastMPdist - distanceMPB);
% figure
% plot(difference)
end