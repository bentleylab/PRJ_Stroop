function chSurrData=buildSurrogate(chData,mode)   
    % Build time/phase shifts surrogate data for each channel
    
    [nS,nT]=size(chData); % nSamples x nTrials
    chSurrData=zeros(nS,nT);
    switch mode
        case 'time'
            for t=1:nT
                timeShift=floor(rand(1)*(nS+1)); % nS+1 because rand works in an open interval
                chSurrData(:,t)=circshift(chData(:,t),timeShift);          
            end
        case 'phase'
            for t=1:nT
                FFT = fft(chData(:,t));
                randPhase=angle(fft(rand(nS,1)));
                randPhase(1)=angle(FFT(1)); %keep positive or negative mean value untuched
                chSurrData(:,t)=ifft(abs(FFT) .* exp(1i*randPhase),'symmetric');
            end
    end
end
