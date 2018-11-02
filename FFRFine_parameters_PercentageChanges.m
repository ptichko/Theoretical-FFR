% %  REFER TO TICHKO AND SKOE (2017) FOR DETAILS ABOUT METHODOLOGY % %
% Paramters to implement a % change in the latencies of the FFR

%Stimulus/ response frequencies
%frequency = [16 17 18 19 21 22 23 25 26 28 29 31 33 35 37 39 41 44 46 49 52 55 58 62 65 69 73 78 82 87 93 98 104 110 117 123 131 139 147 156 165 175 185 196 208 220 233 247 262 277 294 311 330 349 370 392 415 440 466 494 523 554 587 622 659 698 740 831 880];
frequency = [16 17 18 19 21 22 23 25 26 28 29 31 33 35 37 39 41 44 46 49 52 55 58 62 65 69 73 78 82 87 93 98 104 110 117 123 131 139 147 156 165 175 185 196 208 220 233 247 262 277 294 311];

%Latencies: Sub 1
latency = [0 2.6 4.2 7.8 13.6 23.8];

%Scale amplitude
stimamp = frequency./frequency;

%Amplitude values
amp = [1 1 1 2 3 4];

%Drop out frequencies
dropoutfreq = [880 880 880 880 200 100];  %drop out frequenices for each generator. 

%LPF cut-off frequency: Filter applied to the aggregate waveform prior to FFT.
LP = 200; 

%Change latencies by maximum pertange in and in % intervals
maxperc  = 50; %total perc change
percinterv = 10; %interval perc change

%Specify generator(s) to change: combinations of 1 - 6
gen = [4 5];

%Initalize figure parameters and counter
if maxperc == 0 
    colorVec = jet(1);
    counter = 1;
else
colorVec = jet(maxperc/percinterv); %color vector for plotting
counter = maxperc/percinterv;
end

