# Non-Realtime Dual-Channel Audio Processing for Echo Cancellation

This repository contains a MATLAB script (`EchoCanceller.m`) for non-realtime audio processing, focusing on echo cancellation for dual-channel audio files. The script processes WAV files with near-end microphone and far-end speaker signals, applying frequency-domain NLMS-based echo cancellation to produce a clean output with reduced echo.

## Features

- **Echo Cancellation**: Implements frequency-domain NLMS-based echo cancellation with adaptive spectral smoothing to effectively suppress echo in dual-channel audio.
- **Dual-Channel Support**: Processes dual-channel WAV files (near-end microphone and far-end speaker signals).
- **Dynamic Frequency Range**: Adjusts echo frequency range (100 Hz to 4000 Hz) to cover low-frequency echo and speech frequencies.
- **Sampling Rate Support**: Supports 8 kHz and 16 kHz sampling rates with adjustable parameters (`cohRange`, `mufb`).
- **Post-Processing**: Includes residual echo suppression and dual-channel gain normalization to maintain energy levels.
- **Visualization**: Provides detailed plots for analysis, including time-domain waveforms, spectrograms, and energy (dB) plots for both input and output signals.
- **Energy Compensation**: Ensures energy retention through dual-channel gain normalization, preserving natural speech levels.

## Prerequisites

- **MATLAB**: This script requires MATLAB (tested with R2020a and later versions).
- **Audio File**: A dual-channel WAV file with the following specifications:
  - **Channel 1**: Near-end microphone signal (with echo).
  - **Channel 2**: Far-end speaker signal (reference for echo cancellation).
  - **Supported sampling rates**: 8 kHz or 16 kHz.

## Installation

1. Clone the repository to your local machine:
   ```sh
   git clone https://github.com/MetaAcoustix/EchoCanceller.git
   ```
2. Navigate to the repository directory:
   ```sh
   cd EchoCanceller
   ```
3. Ensure MATLAB is installed and added to your system path.

## Usage

1. **Prepare Input File**:
   - Place your dual-channel WAV file (e.g., `original_signal.wav`) in the repository directory.
   - Ensure the file meets the requirements (dual-channel, 8 kHz or 16 kHz sampling rate).

2. **Run the Script**:
   - Open MATLAB and navigate to the repository directory.
   - Open the script `EchoCanceller.m`.
   - Update the `filename` variable in the script to match your input file:
     ```matlab
     filename = 'original_signal.wav';
     ```
   - Run the script:
     ```matlab
     run EchoCanceller.m
     ```

3. **Output**:
   - The processed audio will be saved as `mout_core_original_signal.wav` in the same directory.
   - Two figures will be displayed:
     - **Input Signal Acoustic Plot**: Time-domain waveforms, spectrogram, and energy plot for the input near-end signal.
     - **Output Signal Acoustic Plot**: Time-domain waveforms, spectrogram, and energy plot for the processed near-end signal.

4. **Analyze Results**:
   - Check the output WAV file for echo suppression and speech clarity.
   - Review the spectrogram to confirm echo reduction.
   - Inspect the energy plot to ensure energy levels are preserved.

## Example

### Input File
- **File**: `original_signal.wav`
- **Format**: Dual-channel WAV, 16 kHz sampling rate
- **Channel 1**: Near-end microphone signal (with echo)
- **Channel 2**: Far-end speaker signal

### Run Command
```matlab
run EchoCanceller.m
```

### Output File
- **File**: `mout_core_original_signal.wav`
- **Format**: Dual-channel WAV, 16 kHz sampling rate
- **Channel 1**: Processed near-end signal (echo suppressed)
- **Channel 2**: Original far-end signal (unchanged)

### Visualization
- **Input Spectrogram**: Shows echo and noise in the near-end signal.
- **Output Spectrogram**: Shows reduced echo with preserved speech.

## Testing

To test the script:
1. Use a dual-channel WAV file with known echo (e.g., recorded during a speakerphone call).
2. Run the script as described in the **Usage** section.
3. Verify the output:
   - Listen to `mout_core_original_signal.wav` to confirm echo reduction.
   - Check the spectrogram for echo suppression (reduced energy in echo frequency bands).
   - Ensure speech clarity is maintained and no significant distortion is introduced.

## Contributing

Contributions are welcome! To contribute:
1. Fork the repository.
2. Create a new branch for your feature or bug fix:
   ```sh
   git checkout -b feature/your-feature-name
   ```
3. Make your changes and test thoroughly.
4. Commit your changes with a descriptive message:
   ```sh
   git commit -m "Add feature: your feature description"
   ```
5. Push to your fork:
   ```sh
   git push origin feature/your-feature-name
   ```
6. Open a Pull Request with a detailed description of your changes.

### Potential Improvements
- Add support for additional sampling rates (e.g., 44.1 kHz).
- Integrate noise suppression for specific frequencies (e.g., 440 Hz).
- Enhance distortion handling with low-pass filtering.
- Add non-linear processing (NLP) for residual echo suppression.

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.
