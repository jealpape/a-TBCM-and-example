function SPL = measures_getSPL(signal,fs) % Input signal in Pa.
  windows = 0.05*fs;
  N = floor(length(signal)/windows);
  m = zeros (1,N);
  for count = 1:N
    pos_from = (count-1)*windows+1;
    pos_to = count*windows;
    m(count) = 20*log10(rms(signal(pos_from:pos_to)/2e-5));
  end
  SPL = mean(m);
end
