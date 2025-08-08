function f0 = measures_getf0(signal,fs)
  windows = 0.05*fs;
  N = floor(length(signal)/windows);
  m = zeros (1,N);
  for count = 1:N
    pos_from = (count-1)*windows+1;
    pos_to = count*windows;
    m(count) = subf0(signal(pos_from:pos_to),fs);
  end
  f0 = mean(m);
end

function f0 = subf0(signal,fs)
  lim_inf = ceil(fs/(500));
  lim_sup = floor(fs/(50));
  U = xcov(signal,'unbias');
  U = U(ceil(end/2):end);
  U = (U(lim_inf:lim_sup)-min(U(lim_inf:lim_sup)))/(max(U(lim_inf:lim_sup)) - min(U(lim_inf:lim_sup)));
  [M,P] = findpeaks(U);
  if isempty(P)
    f0 = NaN;
  else
    P = P(find(M>=0.95,1,'first'));
    if isempty(P) 
      f0 = NaN;
    else
      f0 = fs/(P + lim_inf);
    end
  end
end

% 
% function f0 = subf0(signal,fs)
%   N_points = 13;
%   Pxx = abs(fft(signal,2^N_points));
%   Pxx = 20*log10(Pxx(1:floor(end/2))/max(Pxx(1:floor(end/2))));
%   U = xcov(Pxx,'unbias');
%   U = U(ceil(end/2):end);
% 
%   lim_inf = ceil((2^N_points)*50/fs);
%   lim_sup = floor((2^N_points)*500/fs);
% 
%   [M,P] = findpeaks(U(lim_inf:lim_sup));
% 
%   if isempty(P)
%     f0 = 0;
%   else
%     P = P + lim_inf;
%     M = (M-min(M))/(max(M)-min(M));
%     P = P(find(M==1,1,'first'));
%     f0 = P*fs/2^N_points;
%   end
% end