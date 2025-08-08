clc; close all; clear;
%%
addpath("tools/")
addpath("a-TBCM/")

vow = 'AE';
gender = 'male';
resul = ['./QinTA/'];
fs = 44100;

a_PL = [700 1000]; %[300:100:700 750:50:1200 1300:100:1700];
a_Q = 0.5:0.125:1;

ACT = combvec(a_PL,a_Q);
vPL = ACT(1,:);
vQ = ACT(2,:);

a_CT = 0.1:0.2:0.9; 
a_TA = 0.1:0.2:0.9;
a_LCA = [0.4 0.5 0.6];

ACT = combvec(a_CT,a_TA,a_LCA);
vCT = ACT(1,:);
vTA = ACT(2,:);
vLCA = ACT(3,:);

N = length(vCT);
M = length(vPL);

%%
%poolobj = parpool('local',3); 
for n = 1:N
    
    SimResults = struct();
    
    ct = vCT(n);
    ta = vTA(n);
    lca = vLCA(n);
    ia = lca;
    pca = 0;

	H1H2 = zeros(M,1);HRF = zeros(M,1);    MFDR = zeros(M,1);
    ACFL = zeros(M,1);    SQ = zeros(M,1);    OQ = zeros(M,1);
    F0 = zeros(M,1);    SPL = zeros(M,1);    PS = zeros(M,1);
    PC = zeros(M,1);    PL = zeros(M,1);    CPP = zeros(M,1);
    PGO = zeros(M,1);
    
    CT = zeros(M,1);    TA = zeros(M,1);    LCA = zeros(M,1);
    PCA = zeros(M,1);    Q = zeros(M,1);    AsA = zeros(M,1);
    AsP = zeros(M,1);    EJE = zeros(M,1); FOL = zeros(M,1); SIMU = zeros(M,1);
    
    id = 1;
    for a = 1:M                
        q = vQ(a);
        pl = vPL(a);
        act_R = [lca;ia;pca;ct;ta]; act_L = [lca;ia;pca;ct;q*ta];

        [dXR,dXL,Lg,X_L,X_R,a_g,Ut,Pgo,Pout,Psub,Pcol,tm] = Simu_aTBCM(act_R,act_L,pl,vow,gender);

        [asa, asp, eje, h1h2, hrf, mfdr, acfl, f0, sq, oq, cpp, naq] = features(X_L,X_R,Ut,fs);
        spl=measures_getSPL(Pout,fs); ps=mean(Psub); pc=max(Pcol); pgo = Pgo(end);

        % save signals
        SimResults.dXR(:,:,id)=dXR; SimResults.dXL(:,:,id)=dXL; SimResults.Lg(:,:,id)=Lg; SimResults.X_L(:,:,id)=X_L;
        SimResults.X_R(:,:,id)=X_R; SimResults.a_g(:,id)=a_g; SimResults.Ut(:,id)=Ut; SimResults.PGO(:,id)=Pgo;
        SimResults.Pout(:,id)=Pout; SimResults.Psub(:,id)=Psub; SimResults.Pcol(:,id)=Pcol; SimResults.act_R(:,id)=act_R;
        SimResults.act_L(:,id)=act_L; SimResults.PL(id)=pl;

        % save table
        FOL(id) =n; SIMU(id) =id;
        Q(id) = q; PCA(id) = pca; AsA(id) = asa; AsP(id) = asp; EJE(id) = eje;
        PL(id) = pl; CT(id) = ct; TA(id) = ta; LCA(id) = lca; PGO(id) = pgo;
        H1H2(id) = h1h2; HRF(id) = hrf; MFDR(id) = mfdr; ACFL(id) = acfl; CPP(id) = cpp;
        SQ(id) = sq; OQ(id) = oq; F0(id) = f0; SPL(id) = spl; PS(id) = ps; PC(id) = pc;

        id = id + 1;
    end
    
    vec = [FOL SIMU Q, PCA, CT, TA, LCA, PL, H1H2, HRF , MFDR, ACFL , ...
           SQ , OQ ,F0 , CPP, SPL, PS, PC, PGO, AsA, AsP, EJE];
    
    Tab= array2table (vec, 'VariableNames', {'FOL' 'SIMU' 'Q' 'a_PCA' 'a_CT' 'a_TA' 'a_LCA' 'PL' 'H1H2' 'HRF' 'MFDR' 'ACFL' 'SQ' 'OQ' 'F0' 'CPP' 'SPL' 'PS' 'PC' 'PGO' 'AsA' 'AsP' 'EJE'});
    
    parsave([resul 'Table/Tab_' int2str(n)],Tab)
    parsave([resul 'Signal/Data_' int2str(n)],SimResults)
    
end
%delete(poolobj)

