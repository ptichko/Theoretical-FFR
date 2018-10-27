clear all
close all

% %  REFER TO TICHKO AND SKOE (2017) FOR DETAILS ABOUT METHODOLOGY % %

%Stimulus/ response frequencies
frequency = [32 43 51 61 73 87 100 100 155 207 256 293];

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
percinterv = 5; %interval perc change

%Specify generator(s) to change: combinations of 1 - 6
gen = [2 6];

%Initalize figure parameters
colorVec = hsv(maxperc/percinterv); %color vector for plotting

%Figure
figure;
xlim([0 frequency(length(frequency))])
set(gca, 'YTick', []);
ylabel('Relative Amplitude');
xlabel('Frequency (Hz)')
hold on;

%Loop that gradually delays certain generators
for n = 1:(maxperc/percinterv)
    
    %Copy latencies before change by perc
    lantencyperc = latency;
    
    %Multiplier to change latency(ies)
    percmult = 1 + (percinterv * n)/100; 
    
    %Change latency(ies) of generator(s) by percentage
    for g = 1:length(gen)
        lantencyperc(gen(g)) = lantencyperc(gen(g)) * percmult;
    end
    
    %Compute theorectical FFRs
    [frequency, theorNoLP, theor]= generate_TheoFFRTalk(lantencyperc/1000,  dropoutfreq,amp, LP, frequency, stimamp);
    
    %Plot
    plot(frequency,theor, 'ko-', 'LineWidth', 2, 'Color', colorVec(n,:)); 
    
    %Update legend
    legendVec{n}=strcat(num2str(n * percinterv), '%');
   
end

%Update legend
legend(legendVec)

hold off;

