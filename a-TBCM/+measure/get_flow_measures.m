function [H1H2,HRF,MFDR,ac_flow,f0,varargout]=get_flow_measures(ug,fs,varargin)

% try
%Get rid of start and end for signal analysis
n1=round(length(ug)/6);
n2=length(ug)-round(length(ug)/6);
ug=ug(n1:n2);
% ug=ug+1e-4*max(abs(ug)).*randn(1,n2-n1+1);% add some noise for "avoid" zero amplitude in spectra.
dUg=[diff(ug) 0]*fs;

% %Get Spectrum
% nfft=length(ug);
% win=hann(nfft)';
% UG=fft(ug.*win,nfft)/sum(win);
% PSD1=20*log10(abs(UG)/max(abs(UG)));
% PSD1=PSD1(1:round(end/2));
% f=linspace(0,fs/2,length(PSD1));

% figure,
% plot(f,PSD1);
% pause();

% Welch-based periodogram
nfft = max(length(ug),2^10);
Welch_window = length(ug)/2/fs*1000; %ms
win=hann(floor(Welch_window/1000*fs));
[pxx1,f] = pwelch(ug-mean(ug),win,0,nfft,fs,'onesided','ms');
PSD1=10*log10(abs(pxx1)/max(abs(pxx1)));
f=linspace(0,fs/2,length(PSD1));

%Get f0 and period
th_corr=0.5;% between [0 1]
graph='off';%'on';%
displayInfo='off';%'on';%

if nargin<3
    f0= measure.pitchAutoCorr(ug.*hann(length(ug))',fs,th_corr,graph,displayInfo); % used as Simple VAD and f0 estimation
else
    f0=varargin{1};
end



% figure
% pwelch(ug,length(uga0),0,length(uga0))
%
% hold on
% pwelch(uga0,length(uga0),0,length(uga0))
% pwelch(uga1,length(uga1),0,length(uga1))
% hold off
% pause

try
    if isnan(f0)
        % ------------------ Peak estimation -----------------
        [maxtab, mintab, f0]=peakdet(ug, rms(ug),fs);
        % ---------------------------------------------------
        n0=nanmedian(diff(maxtab(:,1)))
        f0=f0(1,1)
    else
        n0=round(fs/f0);
    end
catch
    f0=f(find(PSD1==max(PSD1)))
    n0=round(fs/f0);
end
if isnan(n0)
    disp('f0 is  NaN !!')
    f0=NaN;
end

if not(isnan(f0))
    [b,a]=butter(2,[f0/2 3*f0/2]/(fs/2));
    uga0=filtfilt(b,a,ug);
    a0=rms(uga0);
    
    [b,a]=butter(2,[3*f0/2 5*f0/2]/(fs/2));
    uga1=filtfilt(b,a,ug);
    a1=rms(uga1);
else
    a0=NaN;
    a1=NaN;
end




% % %Simple VAD
% if rms(ug)<30 || f0<80
%     H1H2=inf;
%     HRF=inf;
%     MFDR=inf;
%     ac_flow=inf;
%     f0=inf;
%     varargout{1}=inf;
%     varargout{2}=inf;
%     return
% end

%Loop in time to get temporal measures
m=1;

for k=1:n0:n2-n1-n0
    %Get AC Flow and MFDR
    ac_flow(m) = max(ug(k:k+n0))-min(ug(k:k+n0));
    MFDR(m)    = -min(dUg(k:k+n0));
    if nargout>5
        %Get SQ and OQ
        [dummy,n_max_local]=max(ug(k:k+n0));
        n_max=n_max_local+k;
        n_in=n_max;n_out=n_max;
        threshold=0.5*ac_flow(m);%50% decay to find initial baseline
        while (ug(n_out)>ug(n_max)-threshold)&&(n_out<length(ug))
            n_out=n_out+1;
        end
        while (ug(n_in)>ug(n_max)-threshold) && (n_in>1)
            n_in=n_in-1;
        end
        
        %Linear regression to extend the range
        try
            xL=[n_in,n_max];
            xR=[n_max,n_out];
            yL=[ug(n_in),ug(n_max)];
            yR=[ug(n_max),ug(n_out)];
            pL=polyfit(xL,yL,1);
            pR=polyfit(xR,yR,1);
            mL=pL(1);nL=pL(2);
            mR=pR(1);nR=pR(2);
            n_in_star=round((ug(n_max)-0.95*ac_flow(m)-nL)/mL);
            n_out_star=round((ug(n_max)-0.95*ac_flow(m)-nR)/mR);
            n_in=n_in_star;
            n_out=n_out_star;
            
            %SQ and OQ definitions
            SQ(m) = 100*(n_max-n_in)/(n_out-n_max);
            OQ(m) = 100*(n_out-n_in)/n0;
        catch
            SQ(m) = NaN;
            OQ(m) = NaN;
        end
    end
    m=m+1;
end

cycles=m;
MFDR=nanmean(MFDR);
ac_flow=nanmean(ac_flow);
if nargout>5
    SQ = SQ(3:end-2);
    OQ = OQ(3:end-2);
    SQ_mean=nanmean(SQ);
    OQ_mean=min(nanmean(OQ),100);
    varargout{1}=SQ_mean;
    varargout{2}=OQ_mean;
    varargout{3}=a0;
    varargout{4}=a1;
end
%loop in frequency to get harmonics
m=1;BW=f0*.3;
for k=f0:f0:(fs/2)
    range=find(f>k-BW & f<k+BW);
    harm(m,1)= max(PSD1(range));
    m=m+1;
end

%In case simple VAD did not catch it...
if ~exist('harm','var')
    H1H2=inf;
    HRF=inf;
    MFDR=inf;
    ac_flow=inf;
    f0=inf;
    varargout{1}=inf;
    varargout{2}=inf;
    return
end

H1H2=-harm(2,:);
HRF=10*log10(sum(10.^(harm(2:end,:)/10)));
varargout{5}=rms(ug);

%Get OQ SQ by hand for testing purposes (deactivated)
GUI='off';
if strcmp(GUI,'on')
    figure,
    plot(ug)
    title('ZOOM in to show one full cycle')
    reply = input('Is the zoom adequate? Y/N [Y]: ', 's');
    if isempty(reply), reply = 'Y';end
    title('Select FOUR points: n_{in} n_{max} n_{out} n_{in2}')
    [samples,dummy]=ginput(4);
    n_in=round(samples(1));
    n_max=round(samples(2));
    n_out=round(samples(3));
    n_in2=round(samples(4));
    
    SQ = (n_max-n_in)/(n_out-n_max)*100
    OQ = (n_out-n_in)/(n_in2-n_in)*100
end

Cpp= measure.CPP_HillmanUpdate(ug',f0,fs,hann(length(ug)));%
NAQ=ac_flow*f0/MFDR;

varargout{6}=Cpp;
varargout{7}=NAQ;

% catch error
%   disp(['get_flow_measures_Journal_Z_Scores , ' error.message])
%   H1H2=0;
%   HRF=0;
%   MFDR=0;
%   ac_flow=0;
%   f0=0;
%   varargout={0,0,0};
% end