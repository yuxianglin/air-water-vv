
%------------------------------------------------------
% FUNCTION fft_gauge
% Takes data about wave heights over time, returns FFT decomposition.
% INPUT:
%   - time: array, time each data point was collected.
%   - data_col: height data for one individual gauge over all time.
%   - avg_over: number of values on each side over with to average
% OUTPUT:
%   - f: spectral frequency.
%   - df: increment between two calculated frequencies.
%   - P1: contribution of each frequency, as calculated through FFT.
%------------------------------------------------------

function [f, df, P1] = fft_gauge(time, data_col, avg_over)
    dt = time(2)-time(1);
    tt = time;
    yy = data_col;
    L = length(tt);

    yy = yy - mean(yy);
    
    Y = fft(yy);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = (1/dt)*(0:L/2)/L;
    
    P1( f > 2 ) = [];
    f ( f > 2 ) = [];

    % Smooth spectrum a bit
%     P1 = P1  + [0; 0; P1(1:end-2)] + [0; P1(1:end-1)] + [P1(2:end); 0] + [P1(3:end); 0; 0] ;
%     P1 = P1'/5;
    P1new = P1;
    for I=1:avg_over
        P1new = P1new + [zeros(I,1); P1(1:end-I)] + [P1(I+1:end) ;zeros(I,1)];
    end
    P1 = P1new'/(avg_over*2+1);
    
    df = f(2) - f(1);
end