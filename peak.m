function [ freq, mag ] = peak( x )
%PEAK 
%   Returns the frequency and amplitude of
%   the principal frequency component

    Fs = 8000;
    spect = abs(fft(x)/length(x));
    spect = spect(1:ceil(end/2)); % 0 - 4000 Hz
    [mag ifreq] = max(spect);
    freqs = 0:Fs/length(x):Fs-Fs/length(x);
    freq = freqs(ifreq);
end