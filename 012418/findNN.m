function dist = findNN(x,y)

    %x is the data, y is the query
    n = length(x);
    y = (y-mean(y))./std(y,1);                      %Normalize the query. If you do not want to normalize just comment this line.
    m = length(y);
    x(n+1:2*n) = 0;                                 %Append zeros
    y = y(end:-1:1);                                %Reverse the query
    y(m+1:2*n) = 0;                                 %Append zeros

    %The main trick of getting dot products in O(n log n) time. The algorithm is described in [a].
    X = fft(x);                                     %Change to Frequency domain
    Y = fft(y);                                     %Change to Frequency domain
    Z = X.*Y;                                       %Do the dot product
    z = ifft(Z);                                    %Come back to Time domain

    %compute y stats -- O(n)
    sumy = sum(y);
    sumy2 = sum(y.^2);

    %compute x stats -- O(n)
    cum_sumx = cumsum(x);                           %Cumulative sums of x
    cum_sumx2 = cumsum(x.^2);                       %Cumulative sums of x^2
    sumx2 = cum_sumx2(m:n)-[0;cum_sumx2(1:n-m)];      %Sum of x^2 of every subsequences of length m
    sumx = cum_sumx(m:n)-[0;cum_sumx(1:n-m)];         %Sum of x of every subsequences of length m
    meanx = sumx./m;                                %Mean of every subsequences of length m
    sigmax2 = (sumx2./m)-(meanx.^2);
    sigmax = sqrt(sigmax2);                         %Standard deviaiton of every subsequences of length m


    %computing the distances -- O(n) time. The formula is described in [b].
    dist = (sumx2 - 2*sumx.*meanx + m*(meanx.^2))./sigmax2 - 2*(z(m:n) - sumy.*meanx)./sigmax + sumy2;
    dist = abs(sqrt(dist));

    %If you want Pearson's correlation coefficients instead of Euclidean
    %Distance uncomment the next line. The formula is described in [c].
    %CorrCoef = 1-abs(dist)./(2*m);
end
