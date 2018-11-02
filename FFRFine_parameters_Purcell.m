% % REFER TO PURCELL ET AL (2004) and TICHKO AND SKOE (2017) FOR DETAILS ABOUT METHODOLOGY % %
% Model parameters from Purcell et al. (2004), JASA, to generate EFRs
% 2-generator model with one subcortical G and one cortical G

%Modulation frequencies
frequency = 20:1:100;

%Latencies: Purcell Model
latency = [0 7.3 0 0 29 0];

%Scale amplitude
stimamp = frequency./frequency;

%Amplitude values
% G > 4 will have amplitude ramped down linearly, so place G2 above G4
amp = [0 0.35 0 0 0.85 0];

%Drop out frequencies
dropoutfreq = [110 110 110 110 110 110];  %drop out frequenices for each generator. 

%LPF cut-off frequency: Filter applied to the aggregate waveform prior to FFT.
LP = 110; 


