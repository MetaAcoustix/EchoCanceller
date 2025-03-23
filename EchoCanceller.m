clc
clear

% File name definition
filename = 'original_signal.wav';  
outfilename = ['mout_core_' filename];

% Read audio file and validate
try
    [sig, fs] = audioread(filename);
    if size(sig, 2) < 2
        error('Input audio must be dual-channel');
    end
catch ME
    error('Unable to read audio file: %s', ME.message);
end

% Signal assignment
rrin = sig(:,2);  % far-end spk
ssin = sig(:,1);  % near-end mic

% Basic parameter settings
N = 128;          % Frame length
mult = fs/8000;
if fs == 8000
    cohRange = 2:3;
    mufb = 0.6;
elseif fs == 16000
    cohRange = 2;
    mufb = 0.5;
else
    error('Unsupported sampling rate: %d Hz', fs);
end

% Dynamically adjust echo frequency range
lowFreq = 100;  % Cover more low-frequency echo
highFreq = min(4000, fs/2);  % Cover more speech frequencies
echoBandRange = ceil(lowFreq*2/fs*N):floor(highFreq*2/fs*N);

% Algorithm parameters
AECon = 1;           % AEC switch
FFT_Factor = N;      % FFT scaling factor
suppState = 1;       % Suppression state
len = length(ssin);

% Check input length
if len < N
    error('Input signal length (%d) is less than minimum processing length (%d)', len, N);
end

% Initialize dual-channel output
ercn = zeros(len, 2);
ercn(:,2) = rrin;    % Retain original far-end signal
zm = zeros(N,1);
NN = len;
Nb = floor(NN/N);
xo = zeros(N,1);
do = xo;
eo = xo;

% Tracking variables
hnled = zeros(N+1,1);
weight = zeros(N+1,1);
hnl = zeros(N+1,1);
mbuf = zeros(2*N,1);

% Initialize parameters
hnlLocalMin = 1;
cohxdLocalMin = 1;
ovrd = 2;
ovrdSm = 2;
hnlMin = 1;
hnlMinCtr = 0;
hnlNewMin = 0;

% Window function
wins = [0; sqrt(hanning(2*N-1))];
cohxd = zeros(N+1,1);
Sd = zeros(N+1,1);
Sx = zeros(N+1,1);
Sxd = zeros(N+1,1);

% Processing progress bar
hid = waitbar(0, 'Processing...');
tic

% Main processing loop
for kk = 1:Nb
    pos = N * (kk-1) + 1;
    if pos+N-1 > len
        break;
    end
    
    xk = rrin(pos:pos+N-1);
    dk = ssin(pos:pos+N-1);
    
    xx = [xo; xk];
    xo = xk;
    dd = [do; dk];
    do = dk;
    
    % Adaptive smoothing factor
    signalEnergy = mean(abs(dk).^2);
    gamma = 0.9 + 0.05 * (signalEnergy > 0.01);
    
    tmp = fft(xx .* wins) / FFT_Factor;
    xf = tmp(1:N+1);
    tmp = fft(dd .* wins) / FFT_Factor;
    df = tmp(1:N+1);
    
    Sd = gamma*Sd + (1-gamma)*real(df.*conj(df));
    Sx = gamma*Sx + (1-gamma)*real(xf.*conj(xf));
    Sxd = gamma*Sxd + (1-gamma)*xf.*conj(df);
    
    % Enhanced echo detection
    cohxd = real(Sxd.*conj(Sxd)) ./ (Sx.*Sd + 1e-26);
    cohxd = min(max(cohxd, 0), 1);
    hnled = 1 - cohxd;
    
    % Improved energy detection
    xEnergy = mean(abs(xk).^2);
    dEnergy = mean(abs(dk).^2);
    energyRatio = xEnergy / (dEnergy + 1e-6);  % Far-end to near-end energy ratio
    if xEnergy < 1e-3 || energyRatio < 0.2  % Adjusted threshold
        hnled = ones(N+1,1) * 0.98;  % Light suppression
    end
    
    hnlSortQ = mean(1 - cohxd(echoBandRange));
    [hnlSort2, ~] = sort(hnled(echoBandRange));
    qIdx = floor(0.75 * length(hnlSort2));
    qIdxLow = floor(0.5 * length(hnlSort2));
    hnlPrefAvg = hnlSort2(qIdx);
    hnlPrefAvgLow = hnlSort2(qIdxLow);
    
    if hnlSortQ > 0.9
        suppState = 0;
    elseif hnlSortQ < 0.8
        suppState = 1;
    end
    
    if hnlSortQ < cohxdLocalMin && hnlSortQ < 0.75
        cohxdLocalMin = hnlSortQ;
    end
    
    if cohxdLocalMin == 1
        ovrd = 3;
        hnled = 1 - cohxd;
        hnlPrefAvg = hnlSortQ;
        hnlPrefAvgLow = hnlSortQ;
    end
    
    if suppState == 0
        hnled = 1 - cohxd;
        hnlPrefAvg = hnlSortQ;
        hnlPrefAvgLow = hnlSortQ;
    end
    
    if hnlPrefAvgLow < hnlLocalMin && hnlPrefAvgLow < 0.6
        hnlLocalMin = hnlPrefAvgLow;
        hnlMin = hnlPrefAvgLow;
        hnlNewMin = 1;
        hnlMinCtr = 0;
    end
    
    if hnlNewMin == 1
        hnlMinCtr = hnlMinCtr + 1;
        if hnlMinCtr == 2
            hnlNewMin = 0;
            hnlMinCtr = 0;
            ovrd = max(log(0.01)/(log(hnlMin + 1e-26) + 1e-26), 5);
        end
    end
    
    hnlLocalMin = min(hnlLocalMin + 0.0008/mult, 1);
    cohxdLocalMin = min(cohxdLocalMin + 0.0004/mult, 1);
    
    ovrd = min(ovrd, 5);
    ovrdSm = 0.99*ovrdSm + 0.01*ovrd .* (ovrd < ovrdSm) + ...
             0.9*ovrdSm + 0.1*ovrd .* (ovrd >= ovrdSm);
    ovrdSm = min(ovrdSm, 5);
    
    % Improved gain calculation
    aggrFact = 0.3;
    wCurve = [0; aggrFact*sqrt(linspace(0,1,N))' + 0.1];
    weight = wCurve;
    
    hnled = weight.*min(hnlPrefAvg, hnled) + (1-weight).*hnled;
    od = ovrdSm * (sqrt(linspace(0,1,N+1))' + 1);
    
    % Adaptive gain lower limit
    minGain = 0.8 + 0.1 * (energyRatio > 0.5);  % Dynamically adjust based on energy ratio
    hnl = max(hnled .* (1 + od * 0.5), minGain);
    
    % Enhanced spectral smoothing
    hnl = 0.6 * hnl + 0.2 * [hnl(1); hnl(1:end-1)] + 0.2 * [hnl(2:end); hnl(end)];
    
    if AECon == 1
        df = df .* hnl;
    end
    
    Fmix = df;
    tmp = [Fmix; flipud(conj(Fmix(2:N)))];
    mixw = wins .* real(ifft(tmp)) * FFT_Factor;
    mola = mbuf(end-N+1:end) + mixw(1:N);
    mbuf = mixw;
    ercn(pos:pos+N-1, 1) = mola;
    
    % Debug information
    if mod(kk, floor(Nb/10)) == 0
        fprintf('Block %d: xEnergy=%.2e, dEnergy=%.2e, energyRatio=%.2f, hnl mean=%.4f\n', ...
            kk, xEnergy, dEnergy, energyRatio, mean(hnl));
    end
    
    if mod(kk, floor(Nb/100)) == 0
        str = sprintf('Processing(%ds): %d%%', round(toc), floor(kk/Nb*100));
        waitbar(kk/Nb, hid, str);
    end
end

% Post-processing: Residual echo suppression
for kk = 1:Nb
    pos = N * (kk-1) + 1;
    if pos+N-1 > len
        break;
    end
    
    xk = rrin(pos:pos+N-1);
    dk = ercn(pos:pos+N-1, 1);  % Use near-end signal after first processing
    
    xx = [xo; xk];
    xo = xk;
    dd = [do; dk];
    do = dk;
    
    tmp = fft(xx .* wins) / FFT_Factor;
    xf = tmp(1:N+1);
    tmp = fft(dd .* wins) / FFT_Factor;
    df = tmp(1:N+1);
    
    Sd = gamma*Sd + (1-gamma)*real(df.*conj(df));
    Sx = gamma*Sx + (1-gamma)*real(xf.*conj(xf));
    Sxd = gamma*Sxd + (1-gamma)*xf.*conj(df);
    
    cohxd = real(Sxd.*conj(Sxd)) ./ (Sx.*Sd + 1e-26);
    cohxd = min(max(cohxd, 0), 1);
    hnled = 1 - cohxd;
    
    xEnergy = mean(abs(xk).^2);
    if xEnergy < 1e-3
        hnled = ones(N+1,1) * 0.99;  % Lighter suppression
    end
    
    hnl = max(hnled, 0.95);  % Lightweight suppression
    
    df = df .* hnl;
    
    Fmix = df;
    tmp = [Fmix; flipud(conj(Fmix(2:N)))];
    mixw = wins .* real(ifft(tmp)) * FFT_Factor;
    mola = mbuf(end-N+1:end) + mixw(1:N);
    mbuf = mixw;
    ercn(pos:pos+N-1, 1) = mola;
end

% Post-processing: Dual-channel gain normalization and energy compensation
inputEnergy = mean(abs(ssin).^2);
outputEnergy = mean(abs(ercn(:,1)).^2);
energyCompensation = sqrt(inputEnergy / (outputEnergy + 1e-6)) * 1.1;  % Enhanced compensation
ercn(:,1) = ercn(:,1) * energyCompensation;
ercn(:,2) = ercn(:,2) * max(abs(ssin)) / max(abs(ercn(:,2)) + 1e-6);

% Cleanup and output
close(hid);
try
    audiowrite(outfilename, ercn, fs);
catch ME
    error('Unable to write output file: %s', ME.message);
end

disp('Processing completed');

% Display acoustic plots
t = (0:length(ssin)-1)/fs;

% Input signal plot
figure('Name', 'Input Signal Acoustic Plot', 'NumberTitle', 'off');
subplot(3,2,1);
plot(t, ssin);
title('Input Near-End - Time Domain Waveform');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
subplot(3,2,2);
plot(t, rrin);
title('Input Far-End - Time Domain Waveform');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
subplot(3,2,[3,4]);
spectrogram(ssin, 256, 250, 256, fs, 'yaxis');
title('Input Near-End - Spectrogram');
colorbar;
subplot(3,2,[5,6]);
plot(t, 10*log10(abs(ssin).^2 + 1e-6), 'b');
title('Input Near-End - Energy (dB)');
xlabel('Time (s)');
ylabel('Energy (dB)');
grid on;

% Output signal plot
figure('Name', 'Output Signal Acoustic Plot', 'NumberTitle', 'off');
subplot(3,2,1);
plot(t, ercn(:,1));
title('Output Near-End - Time Domain Waveform');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
subplot(3,2,2);
plot(t, ercn(:,2));
title('Output Far-End - Time Domain Waveform');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
subplot(3,2,[3,4]);
spectrogram(ercn(:,1), 256, 250, 256, fs, 'yaxis');
title('Output Near-End - Spectrogram');
colorbar;
subplot(3,2,[5,6]);
plot(t, 10*log10(abs(ercn(:,1)).^2 + 1e-6), 'r');
title('Output Near-End - Energy (dB)');
xlabel('Time (s)');
ylabel('Energy (dB)');
grid on;
set(gcf, 'Position', [100, 100, 800, 800]);