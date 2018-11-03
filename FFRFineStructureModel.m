clear all
close all

% %  REFER TO TICHKO AND SKOE (2017) FOR DETAILS ABOUT METHODOLOGY % %


% Parameters

% Subject 1 (Tichko & Skoe, 2017)
% parameters_PercentageChanges

% Purcell et al. (2004)
parameters_Purcell

%Change latencies by maximum pertange in and in % intervals
maxperc  = 25; %total perc change
percinterv = 5; %interval perc change

%Specify generator(s) to change: combinations of 1 - 6
gen = [5];

%Initalize figure parameters and counter
if maxperc == 0 
    colorVec = jet(1); %color vector for plotting
    counter = 1;
else
colorVec = jet(maxperc/percinterv + 1); 
counter = maxperc/percinterv + 1;
end

%Figure
figure;
xlim([frequency(1) frequency(end)])
set(gca, 'YTick', []);
ylabel('FFT Amplitude');
xlabel('Frequency (Hz)')
hold on;

%Delays latencies
for n = 1:(counter)
    
    %Copy latencies before change by perc
    lantencyperc = latency;
    
    %Multiplier to change latency(ies)
    percmult = 1 + (percinterv * (n-1))/100; 
    
    %Change latency(ies) of generator(s) by percentage
    for g = 1:length(gen)
        lantencyperc(gen(g)) = lantencyperc(gen(g)) * percmult;
    end
    
    %Compute theorectical FFRs
    [frequency, theorNoLP, theor]= generate_TheoFFR(lantencyperc/1000,  dropoutfreq,amp, LP, frequency, stimamp);
    
    %Plot
    plot(frequency,theor, 'ko-', 'LineWidth', 2, 'Color', colorVec(n,:)); 
    
    %Update legend
    legendVec{n}=strcat(num2str((n -1) * percinterv), '%');
   
end

%Update legend
legend(legendVec)

hold off;

