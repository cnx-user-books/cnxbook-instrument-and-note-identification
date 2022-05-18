%consumes some samples in a sample data structure (See ProjectD) and upsamples them to twice their original sampling frequency
function out = upSampleSamples2(samples, factor)
    out = samples;
    for i = 1:length(samples)
        out(i) = struct('note',samples(i).note,'sample',upSample(samples(i).sample,factor));
    end
   
        
%doubles the sampling rate of a sample
function upSampledSample = upSample(sample, factor)
    upSampledSample = zeros(1,factor*length(sample));
    for k = 1:length(sample)
        upSampledSample((k-1)*factor + 1) = sample(k); %insert the original signal into every other position of our new signal
    end

    upSampledSample = fft(upSampledSample);
    upSampledSample = [upSampledSample(1:ceil(length(upSampledSample)/4)) , zeros(1,floor(length(upSampledSample)/2) - 1), upSampledSample(floor(3*length(upSampledSample)/4):end)]; %cut out the middle of the spectrum
    upSampledSample = real(ifft(upSampledSample))'; %take the inverse fft (have to take real to remove some imaginary part due to rounding errors)
