function [ ratio ] = powerratio( x, cutoff )
%[ ratio ] = POOWERRATIO( x, cutoff )
%   Gives the ratio of power for frequencies
%   above versus below the cutoff frequency
%   Assuming sampling rate of 8000 Hz.
%   Gunnar Atli Sigurdsson, Nox Medical 2012

    Fs = 8000;
    N = length(x);
    df = Fs / N;
    spect = abs(fft(x)/N).^2;
    spect = spect(1:ceil(end/2)); % 0 - 4000 Hz
    icutoff = floor(cutoff/df);
    super = sum(spect(icutoff:end));
    sub = sum(spect(1:icutoff-1));
    ratio = super/sub;
end
