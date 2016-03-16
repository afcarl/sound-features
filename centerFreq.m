function [ freq ] = centerFreq( x )
%[ freq ] = centerFreq( x )
%   returns the frequency that devides
%   the power spectrum in two equal parts
%   Gunnar Atli Sigurdsson, Nox Medical 2012

    Fs = 8000;
    N = length(x);
    df = Fs / N;
    spect = abs(fft(x)/N).^2;
    spect = spect(1:ceil(end/2)); % 0 - 4000 Hz
    total = sum(spect);
    i = find(cumsum(spect) > total/2,1);
    freq = i*df;
end

