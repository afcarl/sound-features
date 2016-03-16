function [ factor ] = crestFactor( m )
%CRESTFACTOR 
%   Calculates the crest factor for 
%   the given time signal m
%   Gunnar Atli Sigurdsson, Nox Medical 2012
    RMS = sqrt(sum(m.^2)/length(m));
    m = sort(abs(m));
    L99 = m(floor(end*99/100));
    factor = L99/RMS;
end