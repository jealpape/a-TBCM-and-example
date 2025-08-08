function F0 = f0(signal,fs)
  lim_inf = ceil(fs/(1000));
  lim_sup = floor(fs/(20));
  U = xcov(signal,'unbias');
  U = U(ceil(end/2):end);
  U = (U(lim_inf:lim_sup)-min(U(lim_inf:lim_sup)))/(max(U(lim_inf:lim_sup)) - min(U(lim_inf:lim_sup)));
  [M,P] = findpeaks(U);
  
  if isempty(P)
    F0 = NaN;
  else
    P = P(find(M>=0.7,1,'first'));
    if isempty(P)
      F0 = NaN;
    else
      F0 = fs/(P + lim_inf);
    end

    NFFT = pow2(nextpow2(length(signal)/4));
    [Pxx,Fxx] = pwelch(signal,NFFT,[],[],fs,'onesided');
    
    if ~isnan(F0)
      H = Pxx(find(Fxx>=F0,1,'first'));
      if (10*log10(max(Pxx)/H) > 80)
        F0 = NaN;
      end
    end
    
%     figure(10011011);
%       hold off;
%       plot(Fxx,10*log10(Pxx),'b','LineWidth',2);
%       hold on;
%       grid on;
%       xlim([0 2e3]);
%       plot([F0 F0],[min(ylim) max(ylim)],'r-','LineWidth',2,'MarkerSize',10);
%       drawnow;    
  end
end
