function [frequency, finalwaveBP, peak, peakBP ] = generate_TheoFFRTalk(audio, varargin);

%   Required input argumenets:
%   audio           Specify audio-stimulus type, 'wav' or 'sinewave'. Use 'wav' and the audiofile name to load
%                   a .wav audio file. Use'sinewave' to generate sinewave
%                   stimuli at user-specified frequencies.

%   Optional input arguments:
%   latency         A vector of latency values to offset the
%                   phase (in milliseconds) of each generator.
%
%   dropoutfreq     A vector of drop-out frequencies for each generator.
%                   Above the drop-out frequency, that generator will no longer contribute
%                   to the final, summed waveform.
%
%   amp             A vector of amplitudes values for each generator used to scale the
%                   amplitude of the auditory stimulus.

%   LP              A scalar for the cut-off frequency for the low-pass
%                   filter.
%
%   frequency       A scalar to generate a sinewave at a specific frequency
%                   in Hz.
%
%   dB_scaled       A vector to 
%  'wav'            Followed by a file name to load
%  'sinewave'       Followed by a scalar of vector of frequencies to
%                   generate sinewave stimuli for

%   Examples:


%% Initialize parameters
audiotype   = 'sinewave';
latency     = [0 0 0 0 0 0];             % 0 ms latency
dropoutfreq = [880 880 880 880 200 100]; % drop-out frequencies
amp         = [1 1 1 2 3 4];             % amplitude values
LP          = 200;                       % low-pass filter cut-off freq
frequency   = [440];                     % 440 Hz default sinewave
dB_scaled   = [];                        % no-scaling


%% Parse inputs

for i = 1:length(varargin)
    if strcmpi(varargin{i},'sinewave') && length(varargin) > i + 5 && ~ischar(varargin{i+1}) ...
            && ~ischar(varargin{i+2}) && ~ischar(varargin{i+3}) && ~ischar(varargin{i+4}) && ~ischar(varargin{i+5})
        audiotype   = varargin{i};
        latency     = varargin{i+1};
        dropoutfreq = varargin{i+2};
        amp         = varargin{i+3};
        LP          = varargin{i+4};
        frequency   = varargin{i+5};
        dB_scaled   = varargin{i+6};
    end
    
    if strcmpi(varargin{i},'wav') && length(varargin) > i + 5 && ~ischar(varargin{i+1}) ...
            && ~ischar(varargin{i+2}) && ~ischar(varargin{i+3}) && ~ischar(varargin{i+4}) && ~ischar(varargin{i+5})
        audiotype   = varargin{i};
        audiofile   = varargin{i+2};
        latency     = varargin{i+3};
        dropoutfreq = varargin{i+4};
        amp         = varargin{i+5};
        LP          = varargin{i+6};
        act_freq    = varargin{i+7};
        dB_scaled   = varargin{i+8};
        
        % Load in user-specific audiofile
        [y Fs] = audioread(audiofil); % signal and sample rate
        t = 0 : 1/Fs : (length(y)/Fs - 1/Fs); % calculate time steps
        
    end
    
end


%% MASTER LOOP  
% For sinewaves, loop through and generate a sinewave at each frequency
% For audiofile, loop through one

for f = 1:length(frequency)
 
%% Generate auditory stimulus

    switch(audiotype)
        
        case('wav')
     
        % Load in user-specific audiofile
        [y Fs] = audioread(audiofile); % signal and sample rate
        t = 0 : 1/Fs : (length(y)/Fs - 1/Fs); % calculate time steps
     
        %Plot audiofile
        figure;
        hold on;
        ylabel('Amplitude');
        xlabel('Time (Seconds)');
        title('Auditory Stimulus');
        plot(t, y)
        hold off;
        
    case('sinewave')
        %% Generate Sinewave
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
    
    % Scale amplitude if there is a drop-out frequency
    % Currently, only works for sine-wave stimuli
    for L = 1:length(latency)
        
        matrix(:,L) = matrix(:,L).* amp(L);
        
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
    
    finalwave = sum(matrix,2);
    
    T = length(finalwave) / Fs;
    t = 0 : 1/Fs : T;
    t(:,end) = [];
    
     %Plot audiofile
     figure;
     hold on;
     ylabel('Amplitude');
     xlabel('Time (Seconds)');
     title('FFR Waveform');
     plot(t, finalwave)
     hold off;
    
    
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
