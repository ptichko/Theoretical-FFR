% %  REFER TO TICHKO AND SKOE (2017) FOR DETAILS ABOUT METHODOLOGY % %
% Paramters to implement a % change in the latencies of the FFR

%Stimulus/ response frequencies
%frequency = [16 17 18 19 21 22 23 25 26 28 29 31 33 35 37 39 41 44 46 49 52 55 58 62 65 69 73 78 82 87 93 98 104 110 117 123 131 139 147 156 165 175 185 196 208 220 233 247 262 277 294 311 330 349 370 392 415 440 466 494 523 554 587 622 659 698 740 831 880];
frequency = 20:1:100;

%Latencies: Purcell Model
latency = [0 7.3 0 0 29 0];

%Scale amplitude
stimamp = frequency./frequency;

%Amplitude values
% G > 4 will have amplitude ramped down linearly
amp = [0 0.35 0 0 0.85 0];

%Drop out frequencies
dropoutfreq = [0 110 0 0 110 0];  %drop out frequenices for each generator. 

%LPF cut-off frequency: Filter applied to the aggregate waveform prior to FFT.
LP = 100; 

%Change latencies by maximum pertange in and in % intervals
maxperc  = 0; %total perc change
percinterv = 0; %interval perc change

%Specify generator(s) to change: combinations of 1 - 6
gen = [];

%Initalize figure parameters and counter
if maxperc == 0 
    colorVec = jet(1); %color vector for plotting
    counter = 1;
else
colorVec = jet(maxperc/percinterv); 
counter = maxperc/percinterv;
end

