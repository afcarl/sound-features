%% absolute fourier spectrum tests

syms k;
int(1/T*cos(2*pi*(k-1)*t/T)-1/T*cos(2*pi*(k+1)*t/T), t, T/2, T)

plot(fft(sin(2*pi*10*(0:0.01:10))))
plot(abs(fft(sin(2*pi*10*(0:0.01:10)))))
plot(abs(fft(abs(sin(2*pi*10*(0:0.01:10))))))
plot(abs(fft(abs(sin(2*pi*10*(0:0.0001:10))))))
plot(abs(fftshift(fft(abs(sin(2*pi*10*(0:0.0001:10)))))))
plot(abs(fftshift(fft(abs(sin(2*pi*10*(0:0.001:10)))))))
plot(abs(sin(2*pi*10*(0:0.001:10))))
plot(abs(sin(2*pi*10*(0:0.001:1))))
plot(abs(sin(2*pi*10*(0:0.001:0.1))))
hold all;
plot(abs(sin(2*pi*10*(0:0.001:0.1))))
clf;



plot(2/pi+0*sin(2*pi*10*(0:0.001:0.1))-4/pi/3*cos(2*pi*20*(0:0.001:0.1))-4/pi/15*cos(2*pi*40*(0:0.001:0.1))-4/pi/35*cos(2*pi*60*(0:0.001:0.1))-4/pi/63*cos(2*pi*80*(0:0.001:0.1))-4/pi/99*cos(2*pi*100*(0:0.001:0.1))-4/pi/143*cos(2*pi*120*(0:0.001:0.1)))
hold all;
plot(abs(sin(2*pi*10*(0:0.001:0.1))))

%%

HDR = sopen('testingBiosig1.edf', 'r', [1 2], 'OVERFLOWDETECTION:OFF'); 

Fs = 8000;
step = 60; %2; %[s]
pos = 1000; %3002;

[s,HDR] = sread(HDR, step, pos);
m = s(:,1)/2^15; %norm amplitude

%%

clf; hold all;
plot(abs(fftshift(fft(m/max(m)))))
XLim = get(gca, 'XLim');
YLim = get(gca, 'YLim');
plot(abs(fftshift(fft(abs(m/max(m))))))
%plot(abs(fftshift(fft((m/max(m)).^2))))
axis([XLim YLim]);
axis([2e5 2.8e5 1e2 1e4])
set(gca, 'YScale', 'log')


%%
sets = 1000;
setSize = 1000;

clf;
for i=1:3
subplot(1, 3, i);
hold all;
switch i
    case 1
        asdf = randn(sets, setSize).*repmat(logspace(-2, 6, sets)', 1, setSize);
        title('Gaussian white noise')
    case 2
        asdf = randn(sets, setSize).*(logspace(-2, 2, sets)'*logspace(-2, 2, sets));
        title('Non-Stationary Gaussian white noise')
    case 3
        asdf = randn(sets, setSize);
        asdf(:, 1) = 1000*ones(setSize,1);
        asdf = asdf.*repmat(logspace(-2, 6, sets)', 1, setSize);
        title('Impulse + Gaussian white noise')
end

x = sum(abs(asdf), 2)/setSize; %mean(abs(x))
x2 = sqrt(sum(asdf.^2, 2)/setSize)/sqrt(setSize); %RMS/sqrt(N)
x3 = sqrt(sum(asdf.^2, 2)/setSize); %RMS
[xs I] = sortrows(x);
x2s = x2(I);
x3s = x3(I);

plot(xs)
plot(x2s)
plot(x3s)
legend('Abs sum', 'Scaled RMS', 'RMS')
set(gca, 'YScale', 'log')
%plotPretty
set(findall(0,'Type','Axes'), 'YGrid', 'off')
end
print -depsc absoluteSumBoundaries

%% gaussian noise

clf; hold all;
asdf = randn(sets, setSize).*repmat(logspace(-2, 6, sets)', 1, setSize);
x = sum(abs(asdf), 2)/setSize;
x2 = sqrt(sum(asdf.^2, 2)/setSize);
plot(x, x2, '.')
C = polyfit(x, x2, 1);
plot(xs, polyval(C, xs))
%plot(xs, polyval([0.04 0], xs))
xlabel('\Sigma |x|/N')
ylabel('RMS')
title('Gaussian Noise')
text(1e6, 2e6, 'Slope: 1.25')
plotPretty
print -depsc gaussNoiseAbsSum
