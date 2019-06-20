clear all
close all

%% AUDIO EXAMPLE

% Parameters
audiofile = 'sinesAfade_9_audend.wav';      % audofile to use as stimulus
latency = [0 2.6 4.2 7.8 13.6 23.8]/1000;   % generator latencies (convert from ms to s) 
amp = [1 1 1 2 3 4];                        % amplitude for each generator
dropoutfreq = [880 880 880 880 200 100];    % drop out frequenices for each generator. ** Doesn't work for audiofiles ** 
LP = 200;                                   % low-pass filter cut-off freq
dB_scale = 1;                               % scale amp of stimuli **** Doesn't work for audiofiles ** 


% Compute theorectical FFRs
[frequency, theorNoLP, theor]= generate_TheoFFR('audio', 'wav', audiofile, ... 
    latency,  dropoutfreq, amp, LP, frequency, dB_scale);

% Plot FFR waveform
figure;
hold on;
ylabel('Amplitude');
xlabel('Time (Seconds)');
title('FFR Waveform');
plot(t, finalwave)
hold off;


%% SINEWAVE EXAMPLE %%

% Parameters
frequency = [30:10:300];                    % stimulus frequencies
latency = [0 2.6 4.2 7.8 13.6 23.8]/1000;   % generator latencies (convert from ms to s) 
amp = [1 1 1 2 3 4];                        % amplitude for each generator
dropoutfreq = [880 880 880 880 200 100];    % drop out frequenices for each generator. 
LP = 200;                                   % low-pass filter cut-off freq
dB_scale = frequency./frequency;            % scale amp of stimuli

% Figure of FFT amplitudes
figure;
xlim([frequency(1) frequency(end)])
set(gca, 'YTick', []);
ylabel('FFT Amplitude');
xlabel('Frequency (Hz)')
hold on;

% Compute theorectical FFRs
[frequency, theorNoLP, theor]= generate_TheoFFR('audio', 'sinewave', latency,  dropoutfreq, amp, LP, frequency, dB_scale);

% Plot FFT amplitude
plot(frequency,theor, 'ko-', 'LineWidth', 2); 

hold off;

