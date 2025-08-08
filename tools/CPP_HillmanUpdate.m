function Cpp=CPP_HillmanUpdate(x,f0,fs,w)

l=length(x);
x=x(:);
if isnan(f0)
    n=round(fs/400):round(fs/100);
else
    n=[(round(fs/f0)-10):(round(fs/f0)+10)];
end

c=rceps(w.*normalize(HP_Filter(x,fs)));
c=c(1:floor(l/2));
dBc=dB20(c);
[maxdBc,maxIndexdBc]=max(dBc(n));


%       plot(abs(c)),
%       hold on,
%       plot(n,(c(n)),'r'),
%       plot(n(1)+maxIndexdBc-1,abs(c(n(1)+maxIndexdBc-1)),'ro')
%       title(num2str(maxdBc))
%       hold off
%       pause

% figure
% plot(dBc)
% pause
try
    [b , ~]=robustfit((n(1):floor(l/2))',dBc(n(1):floor(l/2)),'huber');
%     b
%     plot(dBc)
%     hold on
%     plot(n,dB20(c(n)),'r')
%     hold off
%     line([n(1) 256],[b(2)*n(1)+b(1) b(2)*256+b(1)],'color',[1 0 0])
%     pause
    Cpp=maxdBc-(b(2)*maxIndexdBc+b(1));
catch
    Cpp=NaN;
end
