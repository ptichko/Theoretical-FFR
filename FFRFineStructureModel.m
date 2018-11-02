clear all
close all

% %  REFER TO TICHKO AND SKOE (2017) FOR DETAILS ABOUT METHODOLOGY % %
% Parameter scripts

% Vary latencies by percentage
% FFRFine_parameters_PercentageChanges

% Parameters for Purcell et al. (2004)
FFRFine_parameters_Purcell

%Figure
figure;
xlim([0 frequency(length(frequency))])
set(gca, 'YTick', []);
ylabel('Relative Amplitude');
xlabel('Frequency (Hz)')
hold on;

%Loop that gradually delays certain generators
for n = 1:(counter)
    
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

