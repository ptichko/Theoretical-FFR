function [frequency, peak, peakBP, t, finalwaveBP, xFFT, fft_truncBP] = generate_TheoFFR(varargin)


%   Required input argumenets:
%   audiotype       First agument must be 'wav' or 'sinewave'.
%                   Use 'wav', followed by file name, to load a .wav audio file.
%                   Use'sinewave' to generate sinewave stimuli at user-specified frequencies.
%                   
%   latency         A vector of latency values to offset the
%                   phase (in seconds) of each generator.
%
%   dropoutfreq     A vector of drop-out frequencies for each generator.
%                   Above the drop-out frequency, that generator will no longer contribute
%                   to the final, summed waveform. Works only for sinewave
%                   stimuli.
%
%   amp             A vector of amplitudes values for each generator used to scale the
%                   amplitude of each generator.

%   LP              A scalar for the cut-off frequency for the low-pass
%                   filter that is applied to the final waveform.
%
%   frequency       A scalar to generate a sinewave at a specific frequency
%                   in Hz.
%
%   dB_scaled       A vector to scale the amplitude of each stimulus.


%   Output Arguments:
%   frequency       Frequency of auditory stimulus (useful for sinewave)
%   peak            Spectral amplitude of unfiltered FFR waveform
%   peakBP          Spectral amplitude of LPF FFR waveform
%   t               Time step to plot FFR waveform in time domain
%   finalwaveBP     Final FFR waveform in time domain
%   xFFT            Frequency step to plot FFT of LPF FFR waveform
%   fft_truncBP 	Spectrum (FFT) of the LPF FFR waveform

% Examples:

% Audio:
% audiofile = 'sinesAfade_9_audend.wav';      % audofile to use as stimulus
% frequency = 220;                            % Must pass at least one freq to initiate loop in generate_TheoFFR
% latency = [0 2.6 4.2 7.8 13.6 23.8]/1000;   % generator latencies (convert from ms to s) 
% amp = [1 1 1 2 3 4];                        % amplitude for each generator
% dropoutfreq = [880 880 880 880 200 100];    % drop out frequenices for each generator. ** Doesn't work for audiofiles ** 
% LP = 200;                                   % low-pass filter cut-off freq
% dB_scaled = 1;                               % scale amp of stimuli **** Doesn't work for audiofiles ** 
% 
% [frequency, theorNoLP, theor, t, finalwave]= generate_TheoFFR( 'wav', audiofile,...
% latency,  dropoutfreq, amp, LP, frequency, dB_scaled);
% 
% varargin = {'wav', audiofile, latency,  dropoutfreq, amp, LP, frequency, dB_scaled};

% Sinewave:
% frequency = [40:10:300];
% latency = [0 2.6 4.2 7.8 13.6 23.8]/1000; %convert milliseconds to secs
% amp = [1 1 1 2 3 4];
% dropoutfreq = [880 880 880 880 200 100];  %drop out frequenices for each generator. 
% LP = 200; 
% dB_scaled = ones(1, length(frequency));
% generate_TheoFFR('sinewave', latency, dropoutfreq, amp, LP, frequency, dB_scaled)

% varargin = {'sinewave', latency, dropoutfreq, amp, LP, frequency, dB_scaled};

%% Initialize parameters
% audiotype   = 'sinewave';
% latency     = [0 0 0 0 0 0];             % 0 ms latency
% dropoutfreq = [880 880 880 880 200 100]; % drop-out frequencies
% amp         = [1 1 1 2 3 4];             % amplitude values
% LP          = 200;                       % low-pass filter cut-off freq
% frequency   = [440];                     % 440 Hz default sinewave
% dB_scaled   = [1];                       % no-scaling


%% Parse inputs

for i = 1:length(varargin)
    
    if strcmpi(varargin{i},'sinewave') && length(varargin) > i + 5 && ~ischar(varargin{i+1}) ...
            && ~ischar(varargin{i+2}) && ~ischar(varargin{i+3}) && ~ischar(varargin{i+4}) ... 
            && ~ischar(varargin{i+5}) && ~ischar(varargin{i+6})
       
        audiotype   = varargin{i};
        latency     = varargin{i+1};
        dropoutfreq = varargin{i+2};
        amp         = varargin{i+3};
        LP          = varargin{i+4};
        frequency   = varargin{i+5};
        dB_scaled   = varargin{i+6};
    end
    
    if strcmpi(varargin{i},'wav') && length(varargin) > i + 5 ...
            && ~ischar(varargin{i+2}) && ~ischar(varargin{i+3}) && ~ischar(varargin{i+4}) ...
            && ~ischar(varargin{i+5}) && ~ischar(varargin{i+6}) && ~ischar(varargin{i+7})
        
        audiotype   = varargin{i};
        audiofile   = varargin{i+1};
        latency     = varargin{i+2};
        dropoutfreq = varargin{i+3};
        amp         = varargin{i+4};
        LP          = varargin{i+5};
        frequency   = varargin{i+6};
        dB_scaled   = varargin{i+7};
        
    end
    
end


%% MASTER LOOP  
% For sinewaves, loop through and generate a sinewave at each frequency
% For audiofile, loop through one

for f = 1:length(frequency)
 
%% Generate auditory stimulus

    switch(audiotype)
        
        case('wav')
     
        % User-specified audiofile
        [y Fs] = audioread(audiofile);          % signal and sample rate
        t = 0 : 1/Fs : (length(y)/Fs - 1/Fs);   % calculate time steps
        
        Fn = frequency(f);                      % Set Fn to complete loop 

        % Plot audiofile
%         figure;
%         hold on;
%         ylabel('Amplitude');
%         xlabel('Time (Seconds)');
%         title('Auditory Stimulus');
%         plot(t, y)
%         hold off;
%         
    case('sinewave')
        %% Generate Sinewave
        
        % Parameters for sine-wave generation
        
       % Duration [s]
        T=0.200;

        % Sample rate [Hz] Supported by SoundCard (16000,48000,96000,192000)
        Fs = 48000;
       
        % samples
        N = T*Fs;
        
        % samples vector
        t = 0 : 1/Fs : T;
        
        % Frequency [Hz]
        Fn = frequency(f);
        
        % Signal
        y = sin(Fn*2*pi*t);

        % new step added 11/28 to correct for intensity differences across stimuli
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

    end
    
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
    
    
    % Insert each sine-wave or audiofile into designated column. Off set by latency
    for L = 1:length(latency)
        
        if numpoints(L) ~= 0
            matrix(numpoints(L):numpoints(L) + length(y) - 1,L) = (y);
        else
            matrix(:,L) = (y);
        end
    end
    
    
    % Scale amplitudes of generators 
    for L = 1:length(latency)

        matrix(:,L) = matrix(:,L).* amp(L);
        
        % Further scale amplitude if there is a drop-out frequency
        % Currently, only works for sine-wave stimuli
        if strcmp(audiotype, 'sinewave')

            if Fn>=dropoutfreq(L)  % if above the drop out frequency then remove from the bank of generators.
                matrix(:,L) = matrix(:,L) * 0;
            end


            if Fn<dropoutfreq(L)  % ramps down the frequency 
                if L>4
                    if Fn>=dropoutfreq(L)/2

                        matrix(:,L) = matrix(:,L) * (2 - (2 * (Fn./dropoutfreq(L))));
                        %Ramp amplitude linearly down beginning at (drop-out frequency)/2

                    end
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
    
    %%% EXTRACT SPECTRAL AMPLITUDES
    switch(audiotype)
        
        case('wav')
            
            peak(f) = fft_trunc(round(Fn));         % Spectral amplitude is meaningless for audiofiles.
            peakBP(f) = fft_truncBP(round(Fn));     % Don't use. Just need code to finish loop.
            
        case('sinewave')
    
            peak(f) = fft_trunc(round(Fn+1));
            peakBP(f) = fft_truncBP(round(Fn+1));
    end
    

end
