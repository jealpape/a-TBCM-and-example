function H1H2=h1h2(signal,F0,fs)
% función de h1h2 que usa el tener f0

NFFT = pow2(nextpow2(length(signal)/4));
[Pxx,Fxx] = pwelch(signal,NFFT,[],[],fs,'onesided');

if isnan(F0)
    H1H2 = NaN;
else
    
    % calculando h1-h2
    PSD1=10*log10(abs(Pxx)/max(abs(Pxx)));

    [~,idx0] = sort(abs(F0-Fxx));
    [~,idx1] = sort(abs(2*F0-Fxx));
    H1 = max(PSD1(idx0(1:3)));
    H2 = max(PSD1(idx1(1:3))); 

    H1H2 = H1-H2;   

    end

end

