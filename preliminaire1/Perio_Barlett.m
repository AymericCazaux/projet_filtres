function y = Perio_Barlett(x, NFFT)
    
    %Nombre de segments
    L = floor(length(x)/NFFT);

    %Initialisation
    y = zeros(1, NFFT);

    for i = 1:L
        segment = x((i-1)*NFFT + 1 : i*NFFT);
        fft_seg = fft(segment, NFFT);
        Perio_x = (abs(fft_seg).^2) /NFFT;
        y = y + Perio_x;
    end
    y = y/L;

end