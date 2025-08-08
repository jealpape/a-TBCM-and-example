function maxflow = measures_getMaxflow(signal,fs)
  windows = 0.05*fs;
  N = floor(length(signal)/windows);
  m = zeros (1,N);
  for count = 1:N
    m(count) = max(signal((count-1)*windows+1:count*windows));
  end
  maxflow = mean(m);
end