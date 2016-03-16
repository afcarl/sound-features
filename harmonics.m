function [ result harmfreq ] = harmonics( x )
%[ significance frequency ] = HARMONICS( signal )
%   Returns the average amplitude for the largest 
%   estimated harmonic sequence divided by the
%   average amplitude for all harmonic sequences
%   Gunnar Atli Sigurdsson, Nox Medical 2012

    Fs = 8000; %fixed sampling freq
    N = length(x);
    df = Fs/N; % frequency resolution
    dt = 1/Fs;
    train = 4; %harmonic seq must be more than 3 values;
    w = hamming(N);
    spect = abs(fft(w.*x)/N).^2; %power spectra
    
    % added 15.6.2012, cepstrum approach
    max1stharm = 300; %[Hz]
    min1stharm = 30; %[Hz]
    cepstrum = real(ifft(log(abs(fft(x.*w)))));
    
    %find cepstrum harmonic estimate
    i0 = floor(Fs/max1stharm);
    iend = ceil(Fs/min1stharm);
    [~, temp] = max(cepstrum(i0:iend));
    ti = i0+temp;
    fundEst = Fs/ti;
    
    %filtering:
    %ord = 10;
    %spect = filter(ones(1,ord)/ord,1,spect);
    %spect = circshift(spect, -ord);
    %noise = sqrt(mean(spect.^2));
    
    spect = spect(1:ceil(end/2)); % 0 - 4000 Hz
    
    low = 50; %lowest harmonic possible [Hz]
    i0 = ceil(low/df); %start with >=50 Hz
    iend = floor(N/2/train); %harmonic freq upper bound
    
    % estimate significance from cepstrum estimate
    ti0 = floor(0.9*fundEst/df);
    tiend = floor(1.1*fundEst/df);
    tpower = zeros(tiend-ti0+1,1);
    for i = ti0:tiend
        cnt = floor(N/2/i);
        seq = i:i:cnt*i;
        tpower(i-ti0+1) = sum(spect(seq))/cnt;
    end
    resultEst = max(tpower)/mean(spect);
    
    power = zeros(iend-i0+1,1);
    for i = i0:iend %test frequency [index]
        cnt = floor(iend/i); %should be N/2 instead of iend?
        seq = i:i:cnt*i;
        vals = spect(seq);
        temp = cumsum(vals)./(1:length(vals))';
        
        if length(temp) < train
            temp = 0;
        else
            temp = temp(3:end);
        end
        
        power(i-i0+1) = max(temp);
    end
    
    [maxpwr temp] = max(power);
    maxi = temp+i0-1; %index for max harmonic
    harmfreq = maxi*df;
    avgpwr = mean(power); % avg pwr for harmonic sequences
    
    result = maxpwr / avgpwr;
    
    %debug:
    p = gcf;
    figure(7);
    plot(spect);
    hold on;
    plot(maxi, spect(maxi),'xr');
    plot(fundEst/df, spect(maxi),'or');
    text(maxi, spect(maxi), sprintf('  %g [Hz]', harmfreq));
    text(iend, max(spect)*0.7, sprintf('%g/%g=%g', maxpwr, avgpwr, result));
    text(iend, max(spect)*0.2, sprintf('significEst=%g', resultEst));
    hold off;
    set(gca, 'YScale', 'log')
    figure(p);

    power = sortrows([(i0:iend)' power],[2 1]);
    fprintf('\n%g, %g, %g\n', power(end-2:end,1)*df)
end
