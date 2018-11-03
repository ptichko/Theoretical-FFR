function [frequency, peak, peakBP ] = generate_TheoFFRTalk(latency, dropoutfreq, amp, LP, act_freq, dB_scaled);

% load FreqAmpRev  % feeds the frequencies and corresponding intensities to use in the model
% frequency = round(act_freq); 
% dB_scaled = dB./max(dB);
frequency = act_freq;

for f = 1:length(frequency)
    
    %Parameters for sine-wave generation
    %Duration [s]
    T=0.200;
    
    %Sample rate [Hz] Supported by SoundCard (16000,48000,96000,192000)
    Fs = 48000;
    %samples
    N = T*Fs;
    %samples vector
    t = 0 : 1/Fs : T;
    %Frequency [Hz]
    Fn = frequency(f);
    %Signal
    y = sin(Fn*2*pi*t);
    %% new step added 11/28 to correct for intensity differences across stimuli
    y = y.*dB_scaled(f);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RAMP THE SIGNAL
    ramptime = 5;  % 5 up, 5 down
    signal = y;
    
    ramp = hann(round((ramptime/1000)*Fs)+1);
    
    if rem(ramptime, 2)==0   % if even
        % find the mid-point of the ramp.
        z = ceil(length(ramp)/2);
        % ramps up for 5 ms  then flat until 5 ms
        ramp = [ramp(1:z);ones(length(signal)-size(ramp,1),1);ramp(z+1:length(ramp))];
        
    else
        % find the mid-point of the ramp.
        z = floor(length(ramp)/2);
        % ramps up for 7.5 ms  then flat until last 7.5 ms
        ramp = [ramp(1:z);ones(length(signal)-size(ramp,1),1);ramp(z+1:length(ramp))];
    end
    
    y = ramp'.*signal;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Calculate # of zeros (zero padding) for each latency, then create
    %matrix
    for L = 1:length(latency);
        
        if latency(L) == 0
            numpoints(L)  = 1;
        else
            numpoints(L) = round(latency(L) * Fs);
            
        end
        
    end
    
    matrix = zeros((length(y) + numpoints(L)), length(latency));
    size(matrix);
    %Insert each sine-wave into designated column. Off set by latency
    
    for L = 1:length(latency)
        
        if numpoints(L) ~= 0
            matrix(numpoints(L):numpoints(L) + length(y) - 1,L) = (y);
        else
            matrix(:,L) = (y);
        end
    end
    
    for L = 1:length(latency)
        
        matrix(:,L) = matrix(:,L).* amp(L);
        
        if Fn>=dropoutfreq(L)  % if above the drop out frequency then remove from the bank of generators.
            matrix(:,L) = matrix(:,L) * 0;
        end
        
        
        if Fn<dropoutfreq(L)  % ramps down the frequency 
            if L>4
            if Fn>=dropoutfreq(L)/2
                %matrix(:,L) = matrix(:,L) * (1-(Fn./dropoutfreq(L)));
                %%this will immediately halve the amplitude instead of
                %%linear ramp
                
                matrix(:,L) = matrix(:,L) * (2 - (2 * (Fn./dropoutfreq(L))));
                %Ramp amplitude linearly down beginning at drop-out frequency/2
                
            end
            end
         
        end
        
        
        
    end
    
    finalwave = sum(matrix,2);
    
    T = length(finalwave) / Fs;
    t = 0 : 1/Fs : T;
    t(:,end) = [];
    
    %     if length(frequency)>20;  %figure is not generated if more than 20 frequencies are included.
    %     %plot(t, finalwave);
    %     %figure
    % %     figure(1)
    % %     subplot(length(frequency),1,f);
    % %     plot(t,finalwave)
    % %         axis([0 inf -5 5]) %scale x and y axis
    % %
    % %     ylabel(Fn)
    % %     set(gca, 'FontSize', 10);
    % %     if f == length(frequency);
    % %         xlabel('Time (ms)');
    %     end
    %     end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% WITH LPF, applied to the final aggregate wave
    
    rate = Fs;
    [b,a] = butter(4, 200/(rate/2),'low');
    finalwaveBP = filter(b,a,finalwave); % or use filtfilt
    
    %     rate = Fs;
    %     [b,a] = butter(3, 2/(rate/2),'high');
    %     finalwaveBP = filter(b,a,finalwaveBP); % or use filtfilt
    %
    %
    
    %%% CALCULATE FFT
    ramp = hann(size(finalwave, 1));
    fftsignal = abs(fft (finalwave.*ramp,Fs));  % FFT with zero padding.
    fftsignal = fftsignal*(2./length(finalwave)); %scale it to microvolts
    fft_trunc = fftsignal(1:round(length(fftsignal)/2)+1);  %only plot up to Nyquist.
    xFFT = [0:1:Fs/2]';   % step size is 1
    
    
    
    fftsignalBP = abs(fft (finalwaveBP.*ramp,Fs));  % FFT with zero padding.
    fftsignalBP = fftsignalBP*(2./length(finalwaveBP)); %scale it to microvolts
    fft_truncBP = fftsignalBP(1:round(length(fftsignalBP)/2)+1);  %only plot up to Nyquist.
    xFFT = [0:1:Fs/2]';   % step size is 1
    hold on;
    
    
    peak(f) = fft_trunc(round(Fn+1));
    peakBP(f) = fft_truncBP(round(Fn+1));
    
    
   
    %%
end
