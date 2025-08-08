function x=LP_Filter(x,fs,varargin)
% Low pass filter for voice analisys
% Victor M. Espinoza, vespinoza@uchile.cl, 2015, Santiago, Chile.

if nargin>2
    fc=varargin{1};
    order=varargin{2};
else
    fc=1100;% from Perkell et. al. 1991
    order=10;
end
[z, p, k]=cheby2(order,40,fc/(fs/2),'low');
[sos,g]=zp2sos(z,p,k);
m=size(sos);
for ii=1:m(1)
    x=filtfilt(sos(ii,1:3),sos(ii,4:6),x);
end
x=x*g*g; % Two times g because filtfilt

% disp(['LP Filter -> order = ' num2str(order) ' , fc = ' num2str(fc) ' Hz'])

% % Option 2 - simple Butterworth
% [b a]=butter(8,fc/(fs/2));
% x=filtfilt(b,a,x);


% % Option 2 - FIR
%
% x=[x; zeros(1,length(bWfilter))'];
% x=filter(bWfilter,1,x);
% x(1:length(bWfilter)/2)=[];
% x((length(x)-length(bWfilter)/2+1):length(x))=[];
% x=x';