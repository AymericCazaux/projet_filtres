function y = Perio_Welch(x, NFFT)
    
    %Nombre de segments
    L = floor(length(x)/floor(NFFT/2));
    L = L-1;

    %Initialisation
    y = zeros(1, NFFT);

    for i = 1:L
        segment = x((i-1)*floor(NFFT/2) + 1 : (i+1)*floor(NFFT/2));
        fft_seg = fft(segment, NFFT);
        Perio_x = (abs(fft_seg).^2) /NFFT;
        y = y + Perio_x;
    end
    y = y/L;

end