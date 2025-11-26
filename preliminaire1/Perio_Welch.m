function Pxx = Perio_Welch(x, NFFT)   

    Nw = floor(NFFT/2);     % longueur des segments
    R  = floor(Nw/2);       % recouvrement 50%
    w  = hann(Nw)';         % fenêtre Hann
    K  = floor((length(x)-Nw)/R) + 1;
    
    Pxx = zeros(1, NFFT);
    
    for k = 1:K
        idx = (k-1)*R + (1:Nw);
        segment = x(idx) .* w;
        Xk = fft(segment, NFFT);
        Pxx = Pxx + abs(Xk).^2 / (sum(w.^2)); % normalisation 
    end

    Pxx = Pxx / K;

end