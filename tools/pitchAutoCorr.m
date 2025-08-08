function [pitch, xRms]= pitchAutoCorr(x,fs,th_corr,graph,displayInfo,varargin)
% Pitch estimation using sample autocorrelation function
% For f0, voiced and stable frame estimation
% script : Victor Espinoza, vespinozac@gmail.com
% updated: October 15, 2016

% upsampling to obtain a "more" accurate f0 estimation
timesUpsample=4;
x=resample(x,timesUpsample,1);
fs=fs*timesUpsample;
xRms=rms(x);
x=x-mean(x);
N=length(x);
x=x'.*tukeywin(N,0.1);


% Only if f0 is in range [fmin fmax]
fmin=50;% in Hz
fmax=500;% in Hz
nmax=min(N,floor(fs/fmin));
nmin=floor(fs/fmax);

% Bandpass filter for avoid ripple
filtering='off';% 'off' because 'x' is already filtered
if strcmp(filtering,'on')
    [b, a]=butter(1,[fmin fmax]/(fs/2));
    x2corr=filtfilt(b,a,x);
else
    x2corr=x;
end

%--------------------------------
rxx=xcorr(x2corr,'coef');
rxx=rxx(N:end);% one-side autocorrelation

n=0:(N-1);
thLine=0.8-n/(N-1);
thLine= thLine';
over_thLine=nmin-1+find(rxx(nmin:nmax)>thLine(nmin:nmax));


maximo=find(rxx==max(rxx(nmin:nmax)));% first positive peak after center lag
% minimo=find(rxx==min(rxx(1:nmax)));% first negative peak after center lag
% pause

% x0=HP_Filter(x,fs);
% [b a]=butter(1,1000/(fs/2));
% x=filtfilt(b,a,[diff(x) 0]);
x=normalize(x);

if strcmp(graph,'on')
    figure
    subplot 211
    plot(1:N,rxx)
    line([1 N],[th_corr th_corr],'color',[0 0 0])
    line([nmin nmin],[-1 1],'color',[0.5 0.5 0.5])
    line([nmax nmax],[-1 1],'color',[0.5 0.5 0.5])
    hold on
    plot(1:N,thLine,'m')
    plot(over_thLine,rxx(over_thLine),'.r')
    hold off
    ylim([-1 1])
    xlim([1 N])
    title('Auto-correlation (one side) - pitchAutoCorr function')
    xlabel('sample lag')
    ylabel('normalized value')
    
    subplot 212
    plot(x)
    ylim([-1 1])
    
    disp('Paused...')
    pause
    close
    
end
%
%
% x1=x;%filtfilt([1 -1],1,x);
% x2=x1(floor(N*.1):floor(N*0.5));
% x3=x1(floor(N*0.5):floor(N*0.9));
%
% xmean=[mean(x2) mean(x3)];
% xstd=[std(x2) std(x3)];
%
% xstdmax=max(xstd);
% xstd=xstd/xstdmax;
%
%
% xmean02=xmean(1);
% xstd02=xstd(1);
% xmean03=xmean(2);
% xstd03=xstd(2);
%
%
% xstdTH=0.90;% between [0,1]
% if and(and(not(or(maximo<=nmin,maximo>nmax)),rxx(maximo)>th_corr),sum(xstd)>xstdTH*2)
%     f0=fs/maximo;
%     %   fmin=0.5*fs/minimo;
%     %
%     %   percent=0.2;
%     %   if not(and(f0<(1+percent)*fmin,f0>(1-percent)*fmin)) % to check if both peaks (negatives & positive) are in the same range
%     %     f0=NaN;
%     %   end
%     if strcmp(graph,'on')
%         %     subplot 211-
%         hold on
%         line(maximo*[1 1],[-1 1],'color','m')
%         %     legend('auto-correlation','threshold',...
%         %       'lag for fmin','lag for fmax','lag for f0',...
%         %       'Location','NorthOutside'...
%         %       )
%         hold off
%         hold on
%         [maxtab, mintab]=peakdet(rxx',rms(rxx),fs);
%         plot(mintab(:,1), mintab(:,2), 'g*');
%         plot(maxtab(:,1), maxtab(:,2), 'r*');
%         hold off
%
%         subplot 212
%         h(1)=plot(normalize(x),'color',[1 1 1]*0.5);
%         hold on
%         h(2)=plot(x,'b');
%         hold off
%         if sum(xstd)>xstdTH*4
%             set(h(1),'LineWidth',2)
%             set(h(2),'LineWidth',2)
%         else
%             set(h(1),'LineWidth',1)
%             set(h(2),'LineWidth',1)
%         end
%
%         line([floor(N*.1) floor(N*.5)],[xmean02*0+xstd02 xmean02*0+xstd02],'color',[1 0 0])
%
%         line([floor(N*.5) floor(N*.9)],[xmean03*0+xstd03 xmean03*0+xstd03],'color',[1 0 0])
%         ylim([-1 1])
%         %     text(.5,10,['diff mean = ' num2str(diff(xmean))],'HorizontalAlignment','left')
%         %     diff(xmean)
%         xlabel('samples')
%         ylabel('normalized amplitude')
%         legend('X','dX')
%         disp('Paused...')
%         pause
%         %     close
%     end
% else
%     f0=NaN;
% end


if sum(over_thLine)>1
    maximo=find(rxx==max(rxx(over_thLine)),1);
    f0=fs/maximo;
else
    f0=NaN;
end
if strcmp(displayInfo,'on')
    display(['Fundamental frequency = ' num2str(f0) ' Hz'])
end

pitch=f0;