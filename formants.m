function [ F1, F2 ] = formants( x )
%FORMANTS
%   Calculates the principal formants F1 and F2 of the signal x
%   Gunnar Atli Sigurdsson, Nox Medical 2012
    Fs = 8000;
    Fs2 = 2000; %resample to 4000 Hz;
    [B,A] = cheby1(10, 0.1, Fs2/Fs); %nyquist filter
    x = filter(B,A,x);
    x = resample(x, Fs2, Fs); 
    Fs = Fs2;
    
    w = hamming(length(x));
    params = lpc(w.*x,10); %order ~= Fs [kHz] + 5
    b = 1;
    a = [0 -params(2:end)];
    
    % improved McCandless method (clearer peaks)
    r1 = roots(a);
    r2 = r1(imag(r1) > 0);
    dif = abs(diff(r2));
    [~,i] = min(dif);
    radius1 = abs(r2(i));
    radius2 = abs(r2(i+1));
    ra = roots(a)*2/(radius2+radius1);
    rb = roots(b);
    [bnew, anew] = zp2tf(rb,ra,1);
    
    [H, W] = freqz(b,a); %dont use McCandless improvement

    % find maxima
    y = abs(H);
    dy = [0; y(3:end)-y(1:end-2); 0]; %1st central diff
    zc = [abs(diff(dy > 0)); 0]; %is zero-crossing?
    ddy = [0; y(3:end)-2*y(2:end-1)+y(1:end-2); 0];  %2nd central diff
    maxima = (ddy < 0) & zc; %maxima must be concave
    peaks = find(maxima); %approx extreme value indices
    
    % finding maxima, method 2
    r1 = roots([0 -params(2:end)]);
    r1 = r1(imag(r1) >= 0); %remove duplicates
    sortme = [(1:length(r1))' abs(log10(abs(r1)))];
    sortme = sortrows(sortme);
    r1 = r1(sortme(:,1));
    F1 = angle(r1(1))/2/pi*Fs;
    F2 = angle(r1(2))/2/pi*Fs;
    if F2 < 50
       F2 = angle(r1(3))/2/pi*Fs; 
    end
    
    % finding maxima, method 3
    % trim leading and trailing zeros?
    a = [0 -params(2:end)];
    a2 = a(2:end)/a(2); %trim trailing zero and normalize with leading coeff
    n = length(a2);
    c = diag(ones(1,n-2), -1);
    c(:,end) = fliplr(-a2(2:end)); %companion matrix
    %roots = eig(c);
    
    % use power method to solve for dominant eigenvalue?
    % roots are within unit circle, so dominant root should be
    % closest to the unit circle.
    
    
%     %sort by amplitude
%     sortme = [y(peaks) peaks];
%     sortme = sortrows(sortme);
%     peaks = sortme(:,2);
%     
%     if length(peaks) < 2
%         F1 = 0;

%         F2 = 0;
%     else
%         F = W(peaks(end:-1:end-1))*Fs/2/pi;
%         %F = sort(F); %sort by frequency?
%         F1 = F(1);
%         F2 = F(2);
%     end
    
    p = gcf;
    figure(2); clf;
    freqz(b,a);
    title('Formants, simple LPC estimate')
    %if length(peaks) >= 2
    %    subplot(2,1,1);
    %    text(W(peaks(end))/pi, 20*log10(H(peaks(end))), 'F1')
    %    text(W(peaks(end-1))/pi, 20*log10(H(peaks(end-1))), 'F2')
    %end

    figure(3); clf;
    freqz(bnew,anew);
    title('Formants, imp. McCandless estimate')
    
    figure(4); clf;
    xfft = 20*log10(abs(fft(x))/length(x));
    xfft = xfft(1:ceil(end/2));
    plot(0:Fs/2/length(xfft):Fs/2-Fs/2/length(xfft), xfft)
    hold on; 
    plot(F1, max(xfft), 'x');
    plot(F2, max(xfft), 'o');
    hold off;
    title('Formant estimated sound [dB]'), xlabel('freq')
    figure(p)
end
