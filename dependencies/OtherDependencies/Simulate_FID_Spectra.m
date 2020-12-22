function [FID,Spectrum] = Simulate_FID_Spectra(Chemshift,DeltaFrequency,phase0,AcqDelay,T2,S_0,SNR,dwelltime,vecSize,LarmorFreq,SmoothFIDSpan)
%
% Simulate_FID_Spectra Simulates FIDs.
%
% This function was written by Bernhard Strasser, [month] [year].
%
%
% The function simulates an FID and its spectrum
% 
%
%
% [FID,Spectrum] = Simulate_FID_Spectra(Chemshift,DeltaFrequency,phase0,AcqDelay,T2,S_0,SNR,dwelltime,vecSize,LarmorFreq,SmoothFIDSpan)
%
% Input: 
% -         Chemshift                   ...    Chemical Shift. Unit: [ppm].
% -         DeltaFrequency              ...    If you want to measure not centered around e.g. water. Set to 4.65 for water. Unit: [ppm].
% -         phase0                      ...    Phase of signal. Unit: [degree].
% -         AcqDelay                    ...    The acquisition delay, i.e. the time after which the FID is measured, after excitation. Results in first order phase. Unit: [s].
% -         T2                          ...    The T2 relaxation constant of the signal. Unit: [s].
% -         S_0                         ...    The signal at t=0. Unit: [IU].
% -         SNR                         ...    The simulated SNR of the signal. Unit: [1].
% -         dwelltime                   ...    The dwelltime of the signal. Unit: [s]
% -         vecSize                     ...    The vector size, i.e. the number of measured time points. Unit: [1]
% -         LarmorFreq                  ...    The Larmor frequency, i.e. the frequency of TMS or dss. Unit: [Hz].
% -         SmoothFIDSpan               ...    If the FID shall be smoothed in the time domain to increase SNR. Set 0 for no smoothing. Unit: [1].
%
%
% Output:
% -         FID                         ...     The resulting time domain data.
% -         Spectrum                    ...     The fft of FID.
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None

% Further remarks: 


%% Comments, Definitions, Initializations

if(~exist('SNR','var'))
    SNR = 0;
end



%% Define frequency and time

omega = LarmorFreq * (1 + (Chemshift - DeltaFrequency)/10^6) * 2*pi;

% Define time
t_end = AcqDelay + dwelltime*(vecSize - 1);
t=AcqDelay:dwelltime:t_end;




%% Create Noise

% Define Standard Deviation to get SNR in time domain
if(SNR > 0)
    std_FID = S_0/(2*SNR);       %SNR_spectrum = SNR / sqrt(numel(t));
else
    std_FID = 0;
end

% Create gaußian noise with std std_FID and mean 0
Noise = std_FID*complex(randn([1 vecSize]), randn([1 vecSize])); % Is real and imaginary noise uncorrelated?




%% Create FIDs

FID = [t; S_0 * exp(-(omega - LarmorFreq*2*pi) * 1i*t).*exp(-t/T2)*exp(1i*deg2rad(phase0)) + Noise];



%% Smooth FID


if(exist('SmoothFIDSpan','var') && (SmoothFIDSpan ~= 1))
   FID_smooth = transpose(smooth(FID(2,:),SmoothFIDSpan));
   FID = [t; FID_smooth]; 
end





%% Compute PPM-Scale

Chemshift = compute_chemshift_vector(LarmorFreq,dwelltime,numel(FID(2,:))) - (4.65 - DeltaFrequency);     % Dont ask me why we need to subtract 4.65... Somehow in my code water should be always centered...


%% FFT
Spectrum = [Chemshift; fftshift(fft(FID(2,:))) / sqrt(numel(t))];




