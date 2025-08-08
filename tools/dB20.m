function y=dB20(x)

%y=20*log10(abs(x)+1e-6*rand(1,1));
y=20*log10(abs(x));
%disp('dB20 function noise floor = 1e-6 -> -120 dB')