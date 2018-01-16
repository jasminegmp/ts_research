function [x, y] = ts_generator(Fs, L, noise)
    
    %T = 1/Fs;                   % Sample time
    %T = Fs
    %t = (0:L-1)*T;              % Time vector
    ts = Fs;
    x = 0:ts:L;
    
    offset = 4;
    
    %% Random data with piecewise functions
    for k = 1:length(x)
        if x(k) < L/4
            y(k) = sin(x(k)) + (randn(size(k))/noise) + offset;
        elseif x(k) >= L/4 && x(k) < L*5/8
            y(k) = abs(sin(x(k))) + (randn(size(k))/noise) + offset;
        elseif x(k) > L*.64 && x(k) < L*3/4
            y(k) = abs(sawtooth(x(k))) + (randn(size(k))/noise) + offset; 
        else  
            y(k) = (randn(size(k))/noise) + offset;
        end
    end
end