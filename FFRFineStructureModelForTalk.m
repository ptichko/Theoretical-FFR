clear all
close all

% % % %  REFER TO TICHKO AND SKOE (2017) FOR DETAILS ABOUT METHODOLOGY and
% SETTINGS
frequency = [32 43 51 61 73 87 100 100 155 207 256 293];
% load Subj1
latency = [0 2.6 4.2 7.8 13.6 23.8]; %subj 1

stimamp = frequency./frequency;

amp = [1 1 1 2 3 4];
amp2 = amp*2;
amp3 = amp2;
amp3(1) = 1;
amp4= [1 2 2 4 3 4];
amp5 = [1 1 1 2 7 8];
% amp3 = amp3/1.2;

dropoutfreq = [880 880 880 880 200 100];  %drop out frequenices for each generator. 
LP = 200; % filter applied to the aggregate waveform prior to FFT.

[frequency theorNoLP theor]= generate_TheoFFRTalk(latency/1000,  dropoutfreq,amp, LP, frequency, stimamp);
[frequency theorNoLP2 theor2]= generate_TheoFFRTalk(latency/1000,  dropoutfreq,amp2, LP, frequency, stimamp);
[frequency theorNoLP2 theor3]= generate_TheoFFRTalk(latency/1000,  dropoutfreq,amp3, LP, frequency, stimamp);
[frequency theorNoLP2 theor4]= generate_TheoFFRTalk(latency/1000,  dropoutfreq,amp4, LP, frequency, stimamp);
[frequency theorNoLP2 theor5]= generate_TheoFFRTalk(latency/1000,  dropoutfreq,amp5, LP, frequency, stimamp);


subplot(2,3,3);
plot(frequency,theor, 'ko-', 'LineWidth', 3);  
hold on;
plot(frequency,theor3, 'bo-', 'LineWidth', 3); 
xlim([0 350])
legend('Base Model', 'System Wide Neural Gain but No Cochlear Boost');
set(gca, 'YTick', []);
ylim([0 4])
ylabel('Relative Amplitude');

subplot(2,3,1);
plot(frequency,theor, 'ko-', 'LineWidth', 3);  
hold on;
plot(frequency,theor4, 'go-', 'LineWidth', 3); 
xlim([0 350])
legend('Base Model', 'Subcortical Boost');
set(gca, 'YTick', []);
ylim([0 4])
ylabel('Relative Amplitude');
xlabel('Frequency (Hz)')

subplot(2,3,2);
plot(frequency,theor, 'ko-', 'LineWidth', 3);  
hold on;
plot(frequency,theor5, 'mo-', 'LineWidth', 3); 
xlim([0 350])
ylim([0 4])
legend('Base Model', 'Cortical Boost');
set(gca, 'YTick', []);
ylabel('Relative Amplitude');
xlabel('Frequency (Hz)')
% xlim([34 350])



