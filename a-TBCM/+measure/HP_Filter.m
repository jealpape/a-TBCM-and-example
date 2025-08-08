function x=HP_Filter(x,fs,varargin)


if nargin>2
    fc=varargin{1};
    order=varargin{2};
else
    fc=60;% below f0
    order=6;
end


[b a]=butter(order,fc/(fs/2),'high');
% [z,p,k]=cheby2(10,60,fc/(fs/2),'high');
% [sos,g]=zp2sos(z,p,k);
% m=size(sos);
% for ii=1:m(1)
%   x=filtfilt(sos(ii,1:3),sos(ii,4:6),x);
% end
% x=x*g*g;

% disp(['HP Filter -> order = ' num2str(order) ' , fc = ' num2str(fc) ' Hz'])

x=filtfilt(b,a,x);