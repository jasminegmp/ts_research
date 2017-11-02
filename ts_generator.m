function [t, y] = ts_generator(Fs, L, noise)
    
    T = 1/Fs;                   % Sample time
    t = (0:L-1)*T;              % Time vector

    %% Random data with piecewise functions
    for k = 1 : length(t)
        if k < (length(t)/3)
            y(k) = 2*sin(t(k)) + (randn(size(k))/noise);
        % Removing noisy section for now
        %elseif k >= (length(t)/3) && k < length(t)/2
        %    y(k) = (randn(size(100)));
        elseif k >= (length(t)/3) && k < length(t)/2
            y(k) = 3*abs(sin(t(k))/2) + (randn(size(k))/(noise*2));
        else k < 3*length(t)/4
        %else if k < 3*length(t)/4
            y(k) = 3*sin(3*t(k)) + (randn(size(k))/noise);
        %else
        %    y(k) = (randn(size(100)));    
        end

    end
end