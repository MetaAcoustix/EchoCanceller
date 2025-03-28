% Clear command line and workspace
clc;
clear;

% Real-time audio device setup
fs = 44100;  % Sampling rate (Hz)
N = 128;     % Frame length (samples)

% List available input devices
info = audiodevinfo;
disp('Available Input Devices:');
for i = 1:length(info.input)
    disp([num2str(i) ': ' info.input(i).Name]);
end

% Prompt user to select input device
micDeviceIdx = input('Select the microphone device (enter the number): ');

% Validate input device selection
if micDeviceIdx < 1 || micDeviceIdx > length(info.input)
    disp('Invalid microphone selection. Exiting.');
    return;
end

% Extract and clean device name
micDeviceName = info.input(micDeviceIdx).Name;
micDeviceName = regexprep(micDeviceName, '\s*\(Core Audio\)$', '');

% List available output devices
disp('Available Output Devices:');
for i = 1:length(info.output)
    disp([num2str(i) ': ' info.output(i).Name]);
end

% Prompt user to select output device
outDeviceIdx = input('Select the output device (enter the number, or 0 to disable output): ');

% Setup audio device reader
try
    micReader = audioDeviceReader('SampleRate', fs, 'SamplesPerFrame', N, '‘Device', micDeviceName);
catch ME
    disp('Error setting up audio input device:');
    disp(ME.message);
    return;
end

% Setup audio device writer (if output enabled)
if outDeviceIdx > 0
    if outDeviceIdx > length(info.output)
        disp('Invalid output device selection. Disabling output.');
        outDeviceIdx = 0;
    else
        outDeviceName = info.output(outDeviceIdx).Name;
        outDeviceName = regexprep(outDeviceName, '\s*\(Core Audio\)$', '');
        try
            writer = audioDeviceWriter('SampleRate', fs, 'Device', outDeviceName);
        catch ME
            disp('Error setting up audio output device:');
            disp(ME.message);
            disp('Disabling real-time output...');
            outDeviceIdx = 0;
        end
    end
else
    disp('Real-time output disabled.');
end

% Buffers for recording
rawInputBuffer = [];
inputBuffer = [];
outputBuffer = [];
echoEstimateBuffer = [];

% File names for saving audio
rawInputFile = 'raw_input.wav';
inputFile = 'input_with_echo.wav';
outputFile = 'processed_output.wav';
echoEstimateFile = 'echo_estimate.wav';

% Algorithm parameters
AECon = 1;  % Acoustic Echo Cancellation enabled

% Simulate far-end signal (white noise)
farEndSignal = 0.5 * randn(N * 2000, 1);  % Pre-generate for entire duration

% NLMS adaptive filter parameters
filterLength = 768;  % Kept for echo coverage
mu_base = 0.9;       % Kept for convergence
epsilon = 1e-6;      % Small constant to avoid division by zero
w = zeros(filterLength, 1);  % Filter coefficients
farEndBuffer = zeros(filterLength, 1);  % Far-end signal buffer

% Simulated echo path
echoPath = zeros(100, 1);
echoPath(21) = 0.7;  % Main echo tap
echoPath(51) = 0.3;  % Secondary echo tap

% Processing thresholds
noiseThreshold = 0.015;  % Increased to reduce distortion
dtdThreshold = 0.003;   % Reduced for more sensitive double-talk detection
corrThreshold = 0.05;   % Kept for sensitivity
smoothFactor = 0.95;    % Kept for smooth transitions

% Low-pass filter coefficients
cutoffFreq = 8000;  % Kept to preserve speech
[b_low, a_low] = butter(4, cutoffFreq/(fs/2), 'low');  % 4th-order Butterworth

% Start processing
disp('Starting real-time processing with simulated far-end signal...');
disp('Please speak into the microphone to test the input.');
disp('Recording will stop after 5.8 seconds (2000 frames).');

% Real-time processing loop
frameIdx = 1;
maxFrames = 2000;  % ~5.8 seconds at 44100 Hz with N=128
prevGain = 1;
prevFade = 1;

try
    while frameIdx <= maxFrames
        % Read microphone input
        dk = micReader();
        if size(dk, 2) ~= 1
            error('Microphone input must be single-channel');
        end
        
        % Amplify input (further reduced to minimize over-amplification)
        dk = dk;  % Reduced from 2 to 1.5
        rawEnergy = mean(abs(dk).^2);  % Raw input energy
        
        % Extract far-end signal frame
        startIdx = (frameIdx - 1) * N + 1;
        endIdx = frameIdx * N;
        xk = farEndSignal(startIdx:endIdx);
        
        % Simulate echo
        echoSignal = filter(echoPath, 1, xk);
        dk_with_echo = dk + echoSignal;
        farEndEnergy = mean(abs(xk).^2);  % Far-end signal energy
        
        % NLMS adaptive filter with optimized double-talk detection
        mola = zeros(N, 1);
        echoEstFrame = zeros(N, 1);
        if AECon
            % Pre-compute correlation for the entire frame
            corrCoef = abs(corr(dk_with_echo, xk));
            powerRatio = rawEnergy / (farEndEnergy + epsilon);
            % Dynamic double-talk threshold
            dynamicDtdThreshold = dtdThreshold * max(0.5, 1 - rawEnergy / 0.005);
            dynamicCorrThreshold = corrThreshold * max(0.5, 1 - rawEnergy / 0.005);
            isDoubleTalk = (powerRatio > dynamicDtdThreshold) && (corrCoef < dynamicCorrThreshold);
            
            % Dynamic step size
            if isDoubleTalk
                mu = mu_base * 0.5;  % Increased from 0.2 to 0.3 for better adaptation
            else
                mu = mu_base * 1.0;  % Non-double-talk step size
            end
            
            % Update filter for the entire frame
            for n = 1:N
                farEndBuffer = [xk(n); farEndBuffer(1:end-1)];  % Update buffer
                echoEstimate = w' * farEndBuffer;
                echoEstFrame(n) = echoEstimate;
                mola(n) = dk_with_echo(n) - echoEstimate;
                
                % NLMS update
                normFactor = farEndBuffer' * farEndBuffer + epsilon;
                adaptiveMu = mu / normFactor;
                w = w + adaptiveMu * mola(n) * farEndBuffer;
            end
        else
            mola = dk_with_echo;
        end
        
        % Noise gate with adjusted threshold
        fadeFactor = min(1, abs(mola) / noiseThreshold).^0.5;
        fadeFactor = smoothFactor * prevFade + (1 - smoothFactor) * fadeFactor;
        prevFade = fadeFactor;
        mola = mola .* fadeFactor;
        
        % Apply low-pass filter
        mola = filter(b_low, a_low, mola);
        
        % Dynamic gain control with adaptive upper limit
        signalPower = mean(abs(mola).^2);
        maxGain = min(0.85, 1.0 - 0.7 * rawEnergy / 0.1);  % Adaptive max gain based on input energy
        targetGain = min(maxGain, max(0.7, 0.9 * rawEnergy / (signalPower + 1e-6)));
        gain = 0.98 * prevGain + 0.02 * targetGain;  % Smoother transition
        prevGain = gain;
        mola = mola * gain;
        
        % Output audio (if enabled)
        if outDeviceIdx > 0
            writer(mola);
        end
        
        % Store in buffers
        rawInputBuffer = [rawInputBuffer; dk];
        inputBuffer = [inputBuffer; dk_with_echo];
        outputBuffer = [outputBuffer; mola];
        echoEstimateBuffer = [echoEstimateBuffer; echoEstFrame];
        
        % Debug output every 100 frames
        if mod(frameIdx, 100) == 0
            echoEstEnergy = mean(abs(echoEstFrame).^2);  % Calculate echo estimate energy
            erle = 10 * log10(mean(abs(dk_with_echo).^2) / (mean(abs(mola).^2) + 1e-6));  % ERLE
            residualNoiseEnergy = mean(abs(mola - dk).^2);  % Residual noise energy
            echoCancelError = mean(abs(mola).^2);  % Echo cancellation error energy
            disp(['Frame ' num2str(frameIdx)]);
            disp(['Raw input energy: ' num2str(rawEnergy)]);
            disp(['Far-end energy: ' num2str(farEndEnergy)]);
            disp(['Echo signal energy: ' num2str(mean(abs(echoSignal).^2))]);
            disp(['Echo estimate energy: ' num2str(echoEstEnergy)]);
            disp(['Output energy: ' num2str(signalPower)]);
            disp(['Residual noise energy: ' num2str(residualNoiseEnergy)]);
            disp(['Echo cancellation error: ' num2str(echoCancelError)]);
            disp(['Double-talk detected: ' num2str(isDoubleTalk)]);
            disp(['Power ratio: ' num2str(powerRatio)]);
            disp(['Dynamic DTD threshold: ' num2str(dynamicDtdThreshold)]);
            disp(['Dynamic Corr threshold: ' num2str(dynamicCorrThreshold)]);
            disp(['Correlation coefficient: ' num2str(corrCoef)]);
            disp(['ERLE (dB): ' num2str(erle)]);
            saveAudioFiles(rawInputBuffer, inputBuffer, outputBuffer, echoEstimateBuffer, ...
                rawInputFile, inputFile, outputFile, echoEstimateFile, fs);
            disp(['Saved intermediate files at frame ' num2str(frameIdx)]);
        end
        
        frameIdx = frameIdx + 1;
    end
catch ME
    disp('Processing stopped.');
    disp(ME.message);
    rethrow(ME);
end

% Save final audio files
saveAudioFiles(rawInputBuffer, inputBuffer, outputBuffer, echoEstimateBuffer, ...
    rawInputFile, inputFile, outputFile, echoEstimateFile, fs);
disp('Final files saved.');

% Release audio devices
release(micReader);
if outDeviceIdx > 0
    release(writer);
end

% Helper function to save audio files
function saveAudioFiles(rawInput, input, output, echoEst, rawFile, inputFile, outputFile, echoFile, fs)
    if ~isempty(rawInput)
        rawInputNorm = rawInput / (max(abs(rawInput)) + 1e-6);
        audiowrite(rawFile, rawInputNorm, fs);
        disp(['Raw input saved to ' rawFile]);
    end
    if ~isempty(input)
        inputNorm = input / (max(abs(input)) + 1e-6);
        audiowrite(inputFile, inputNorm, fs);
        disp(['Input (with echo) saved to ' inputFile]);
    end
    if ~isempty(output)
        outputNorm = output / (max(abs(output)) + 1e-6);
        audiowrite(outputFile, outputNorm, fs);
        disp(['Processed output saved to ' outputFile]);
    end
    if ~isempty(echoEst)
        echoEstNorm = echoEst / (max(abs(echoEst)) + 1e-6);
        audiowrite(echoFile, echoEstNorm, fs);
        disp(['Echo estimate saved to ' echoFile]);
    end
end