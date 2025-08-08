function  [AsA, AsP, eje, H1H2, HRF, MFDR, ACFL, f0, SQ, OQ, CPP, NAQ]= features(X_L,X_R,Ut,fs)

[H1H2, HRF, MFDR, ACFL, f0, SQ, ~, ~, ~, ~, CPP, NAQ]=get_flow_measures_Journal_Z_Scores(Ut',fs);


% Calculo de Asimetr�as
KR = min(X_R(:,1),X_R(:,2));
KL = min(X_L(:,1),X_L(:,2));

vec = 1*((KR+KL)<0);
[~,loc] = findpeaks(vec);
[~,loc2] = findpeaks(-vec);

eje = mean(KR(loc));

AP = [];
AA = [];
OQ = [];

m =1;
for i = 1:(numel(loc)-1)
 %%
    
    t1 = loc(i);
    t2 = min(loc2((loc2>t1)));
    t3 = loc(i+1);
    
    KR2 = KR(loc(i):loc(i+1));
    KL2 = KL(loc(i):loc(i+1));

    [DR,TR] = max(KR2);
    [DL,TL] = max(KL2);

    TO = sum(1*((KR2+KL2)>0));
    % aMPLITUDES
    AR = DR-eje;
    AL = DL+eje;

    % asimetrias
    %Asimetría de Fase
    AP(m) = 100*(TR-TL)/TO;
    % Asimetría de Amplitud
    AA(m) = 100*(AL-AR)/(AL+AR);
    % OQ
    OQ(m) = 100*(t3-t2)/(t3-t1);
    
    m = m+1;
end



AsP = mean(AP);
AsA = mean(AA);
OQ = mean(OQ);

end