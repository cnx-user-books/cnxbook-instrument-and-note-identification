%  E, F, F#, G, G#, A, A#, B, C, C#, D, D#, E, F, F#, G (lowest to highest frequency)


%A, A#, B, C, C#, D, D#, E, F, F#, G, G#
%1, 2,  3, 4, 5,  6, 7,  8, 9, 10, 11,12

%6, 7,  8, 9, 10, 11,12, 1, 2,  3, 4, 5
 
%samples = struct('note',{'E0', 'F0', 'Fsharp0', 'G0', 'Gsharp0', 'A0', 'Asharp0', 'B0', 'C0', 'Csharp0', 'D0', 'Dsharp0', 'E1', 'F1', 'Fsharp1', 'G1', 'Gsharp1', 'A1', 'Asharp1', 'B1', 'C1', 'Csharp1', 'D1', 'Dsharp1E2', 'F2', 'Fsharp2', 'G2', 'Gsharp2', 'A2', 'Asharp2', 'B2', 'C2', 'Csharp2', 'D2', 'Dsharp2'}, 'sample',{E0, F0, Fsharp0, G0, Gsharp0, A0, Asharp0, B0, C0, Csharp0, D0, Dsharp0, E1, F1, Fsharp1, G1, Gsharp1, A1, Asharp1, B1, C1, Csharp1, D1, Dsharp1E2, F2, Fsharp2, G2, Gsharp2, A2, Asharp2, B2, C2, Csharp2, D2, Dsharp2});

%01:49| Garion176: Just for your reference...
%01:50| Garion176: Chalameau is F0 through E1, Throat is G1 through Asharp1 and Clarion is B1 through C2
%01:50| Garion176: It's missing E0 and F1
%01:52| Garion176: Heh, and Fsharp1, although I have that.  I just ignored it for some reason


% 18:56| Garion176: It's on the clarinet scale... E is the lowest note on a clarinet, so the numbers relate to the octave with base note E.  I never played the Csharp2 or higher notes; I can if you want me to
% 18:56| Garion176: The 33rd note, starting from E, is C
% 18:57| TigerStorm023: its cool
% 18:57| TigerStorm023: so, what order do put them in if i wnat lowest to highest frequency
% 18:57| TigerStorm023: E0... F0..... etc then wrap around to A0?
% 18:58| Garion176: E0 is lowest frequency... 440 is the frequency of B0 (clarinet pitch is not concert pitch)
% 18:58| TigerStorm023: interesting

% 19:06| Garion176: concert + 2 = clarinet
%   A on concert = B on clarinet
% 19:06| TigerStorm023: ahh

%concert pitch note structure for clarinet:
% samples = struct('note',{...
%       6, 7,  8, 9, 10, 11,12, 1, 2,  3, 4, 5,...
%       6, 7,  8, 9, 10, 11,12, 1, 2,  3, 4, 5, ...
%       6, 7,  8, 9, 10, 11,12, 1, 2, ...
%   }, 'sample',{...
%           E0, F0, Fsharp0, G0, Gsharp0, A0, Asharp0, B0, C0, Csharp0, D0, Dsharp0,...
%           E1, F1, Fsharp1, G1, Gsharp1, A1, Asharp1, B1, C1, Csharp1, D1, Dsharp1,...
%           E2, F2, Fsharp2, G2, Gsharp2, A2, Asharp2, B2, C2});

%legend('E0', 'F0', 'Fsharp0', 'G0', 'Gsharp0', 'A0', 'Asharp0', 'B0', 'C0', 'Csharp0', 'D0', 'Dsharp0', 'E1', 'F1', 'Fsharp1', 'G1', 'Gsharp1', 'A1', 'Asharp1', 'B1', 'C1', 'Csharp1', 'D1', 'Dsharp1', 'E2', 'F2', 'Fsharp2', 'G2', 'Gsharp2', 'A2', 'Asharp2', 'B2', 'C2', 'Csharp2', 'D2', 'Dsharp2');


%Description:
%This program consumes a recording of a song and some data about the song,
%and samples of the insturments being used and produces data indicating
%when, how long, and how strongly each not is being played.
%Arguments:
%signal - signal to process
%chunkSize - size of pieces to break the signal up into before processing -- try to make this as large as possible without using a lot of virual memory
%windowSize - size of windows to break the chunks into to check for notes inside -- much smaller than 4096 and the windowing will distort the spectrum, so keep it at a nice 8192 if possible
%samplingRate - the samplign rate used to record the samples and the signal
%maxNotesPerWindow - the maximum number of notes to check for inside a window of the song-- the system will look at the most likely notes first, and its usually right, but the more note you have it look for, the more accurate its recognition will be
%samples - a special structure for holding the note sample (see above)
%Example:
%out = ProjectD(signal2, 8192 * 16, 8192, 44100, 2, samples2); % runs the analysis
%postProcessor(out,0) % produces meaningful graphs of the output


function out = ProjectD(signal, chunkSize, windowSize, samplingRate, maxNotesPerWindow, samples)
    size = length(signal);    % the number of samples in the song
    windows = nonOverlappingWindowUp(signal, chunkSize);  % cuts the song into several smaller, more manageable, segments.
    numChunks = length(windows(1,:));  % the number of segments the song was cut into
    
    %analyzes each segment, and compiles the final analysis
    one = GetNotes(windows(:,1), windowSize, samplingRate, maxNotesPerWindow, samples);
    results = zeros(length(one(:,1)), length(one(1,:)) , numChunks);
    for i = 1:numChunks
        if(i > 1)
            results(:,:,i) = (GetNotes(windows(:,i), windowSize, samplingRate, maxNotesPerWindow, samples)); %the command to analyze this segment
        else
            results(:,:,1) = one;
        end
            
        PercentComplete = 100 * (i-1)/numChunks;
    end

   %append each segment's analysis together
    resultChunkLength = length(results(:,1,1));
    out = zeros(resultChunkLength * numChunks,length(results(1,:,1)));
    for i = 1:numChunks
        out((i-1)*resultChunkLength + 1: (i)*resultChunkLength,:) = results(:,:,i);
    end
    

%This function uses a heuristic technique to rank each note (A...G)'s likelihood of
%being played inside a given window of the song. The heuristic used is to look at
%the spectrum of the (hanning) window of data; at each harmonic of each note, the 
%height of the spectrum is recorded. A score based on a formula (below) is used
%to rank each note's likleihood. Finally, a matrix containing the ranked order
%of each note's likelihood inside each window is produced.
function out = GetNotes(signal, windowSize, samplingRate, maxNotesPerWindow, samples) % samples = struct('note', {1, 2, 3, etc}, 'sample', {sample1, sample2, etc})

    windowMatrix = windowUp(signal, windowSize);    %compute windows of the data
    numWindows = length(windowMatrix(1,:));
    windowSize = length(windowMatrix(:,1));
    
    hanningWindowMatrix = zeros(size(windowMatrix)); %make these windows into hanning windows
    for i = 1:numWindows
        hanningWindowMatrix(:,i) = windowMatrix(:,i) .* HANN(windowSize);
    end
    
    stFFTSignalMatrix = fft(hanningWindowMatrix);    %compute ffts of these windows
    
    startingAFreq = 110;
    numOctaves = floor((log(samplingRate / 2 / startingAFreq) / log(2^(1/12))) / 12); %find out what notes are playing
    stNoteScoreMatrix = zeros(12,numWindows);
    for window = 1:numWindows
        noteScore = zeros(12,1);
        for octave = 1:numOctaves
            for note = 1:12
                freq = noteToFrequency(note,startingAFreq * octave);
                digitalCenter = round(AnalogToDigitalFrequency(noteToFrequency(note,startingAFreq * octave), samplingRate)  * windowSize + 1); % the frequency of this harmonic to check
               
                lowOffset = round( (digitalCenter - AnalogToDigitalFrequency(noteToFrequency(note-1,startingAFreq * octave), samplingRate) * windowSize) / 2);
                lowBound = digitalCenter - lowOffset;
                highBound = digitalCenter + lowOffset;
                
                noteSize = highBound-lowBound + 1; %the size of the range of frequencies to check
                
                noteScore(note) = noteScore(note) + (norm(stFFTSignalMatrix(lowBound:highBound,window)))^2 / noteSize; %the "special formula" Basically, the idea is to find how high the spectrum is, divided by the area we're checking (because, higher harmonics take up more spectral space, we don't want to give them a higher weight than the lower frequencies)
            end
        end
        stNoteScoreMatrix(:,window) = noteScore;
    end
    
    [sortedNoteScores,sortedNotes] = sort(stNoteScoreMatrix);    %produce the rankings by sorting the rating matrix produced above

    %Next, we want to convolve the (boxcar) windows of the sample with our samples of each note being played, and see which sample matches best. In order to avoid false detections of samples of the same note at different octaves, we use a special formula to analyze the output of the matched filter.
    paddedWindowMatrix = zeroPad(windowMatrix);
    numSamples = length(samples);
    sampleScores = zeros(numWindows, numSamples) - Inf;
    for window = 1:numWindows
        for note = 1:maxNotesPerWindow %figure out how much the top possible notes actually match the waveform of the window (checked with convolution... using fft's of course)
            for sample = 1:numSamples
                if(samples(sample).note == sortedNotes(12 - note + 1,window))
                      sampleScores(window, sample) = log(norm(exp(abs(ConvoMachine(windowMatrix(:,window), samples(sample).sample))) ));	%do the convolution of this sample, and then rank it
			    %this ranking formula was chosen because it emphassises peaks in the output of the matched filter,
			    %while diminishing the effects of partial matches (especially those from higher octave notes).
			    %It was chosen based on an educated trial and error session, to find the best preforming formula,
			    %so there isn't a derivation or mathematical proof that this is the best formula,
			    %it just seems to be the best one we've tried.
               end
            end
        end
    end
    
    out = sampleScores;
            
%This is the convolver. Basically, it consumes two signals, zero pads them so that their lengths are equal,
%and then doubles their size (with more zeros), and then takes the ffts, and returns the ifft of the output.        
function out = ConvoMachine(signalA, signalB)
    lengthA = length(signalA);
    lengthB = length(signalB);
    desiredLength = max(lengthA,lengthB);
    
    if(lengthA < desiredLength)
        signalA = [signalA; zeros(desiredLength - lengthA,1)];
    end
    
    if(lengthB < desiredLength)
        signalB = [signalB; zeros(desiredLength - lengthB,1)];
    end
    
    signalA = [signalA; zeros(desiredLength,1)];
    signalB = [signalB; zeros(desiredLength,1)];
    %ok, so now the signals are both long enough that multiplication in the frequency domain will be normal (non circular) convolution in the time domain
    
    fftA = fft(signalA);    %compute ffts of these windows
    fftB = fft(signalB);    %compute ffts of each note
    
    fftA(1) = 0; %oh, and zero out the DC offset of each signal to avoid that type of confusion
    fftB(1) = 0;
 
    out =  ifft(fftA .* fftB); %multiply and inverse transform out to get the convolution
    
% converts a digitalFrequency (fft index / fft size) and a sampling rate (Hz) into the analog frequency it cooresponds to. digitalFreq must be less than or equal to 1/2
function out = digitalToAnalogFrequency(digitalFreq, samplingRate) 
    out = digitalFreq * samplingRate;

 % converts an analog frequency into its digital frequency (freq/samplingRate)
function out = AnalogToDigitalFrequency(freq, samplingRate)
    out = (freq / samplingRate);

%outputs a number between 1 and 12 corresponding to that frequency's note (0 is A, 1 is A#, 2 is B, 3 is C, etc)
function out = frequencyToNote(frequency) 
    out = mod(round(log(frequency/220) / log(2^(1/12))),12) + 1;

%converts a note (1 to 12) to a frequency with freqA as the frequency of A
function out = noteToFrequency(note, freqA) 
    out = freqA * 2^((note-1)/12);       

%takes windows of the sample at every half-time of the windowSize; if the sample does not devide evenly, the last sample is zero padded
function out = windowUp(sample, windowSize) 
    sampleLength = length(sample);
    window = boxcar(windowSize);
    
    numWindows = floor(sampleLength / windowSize) + 1;
    sample = [sample ; zeros(windowSize,1)]; % pad so last window fits in there.
    
    out = zeros(windowSize, numWindows);
    for i = 1:numWindows
        position = (i-1) * floor(windowSize/2) + 1;
        chunk = sample(position:position + windowSize - 1);
        out(:,i) = window .* chunk;
    end

%breaks the sample into non-overlapping windows of the given size; if the sample does not devide evenly, the last sample is zero padded
function out = nonOverlappingWindowUp(sample, windowSize) %takes windows of the sample at every half-time of the windowSize, zeros are added onto the last window if the sample isn't divided evenly into windows.
    sampleLength = length(sample);
    window = boxcar(windowSize);
    
    numWindows = floor(sampleLength / windowSize) + 1;
    
    sample = [sample ; zeros(windowSize,1)]; % pad so last window fits in there.
    
    out = zeros(windowSize, numWindows);
    for i = 1:numWindows
        position = (i-1) * windowSize + 1;
        chunk = sample(position:position + windowSize - 1);
        out(:,i) = window .* chunk;
    end    

%takes a matrix of windowed samples and zero-pads them to twice their length
function out = zeroPad(windowMatrix)
    windowSize = length(windowMatrix(:,1));
    numWindows = length(windowMatrix(1,:));
    out = zeros(windowSize  * 2, numWindows);
    for i = 1:length(windowMatrix(1,:))
        out(:,i) = [windowMatrix(:,i); zeros(windowSize,1)];
    end

    
    
    
    
    

            
