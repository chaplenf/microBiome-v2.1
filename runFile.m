% clear workspace

clear;
    
%  Initialize COBRA Toolbox
%    initCobraToolbox;
    
global community bucket

% set the variable "model" to a COBRA Toolbox model      
     Neuro=readCbModel('iFC578.xml');
     model1 = Neuro;
     model1.description = 'Nitrosomonas europaea';
     Nwino=readCbModel('iFC579.xml');  
     model2 = Nwino;
     model2.description = 'Nitrobacter winogradskyi';

% Initiatize variable community
%   community:  description of the microbial community and its environment
%       community.species:      species in the community
%       community.mets:         identifiers for substrates and products in bioreactor
%       community.Xfeed:        biomass in the feed (gDCW/L)
%       community.Lfeed:        metabolites in the feed (mmol/L)
%       community.LNMfeed:      metabolites in the feed (mmol/L)
%       community.Gfeed         partial pressure of species in the gas feed (atm) 
%       community.Fin:          liquid flowrate in (L/hr)
%       community.Fout:         liquid flowrate out (L/hr)
%       community.Vin:          gas flowrate in (L/hr)
%       community.Vout:         gas flowrate out (L/hr)

    community.species = {model1 model2};
    community.mets = {'EX_cpd00075_e','EX_cpd00209_e','EX_cpd00418_e','EX_cpd00007_e','EX_cpd00013_e','EX_cpd00047_e','EX_cpd00659_e','EX_cpd00528_e','EX_cpd00011_e','EX_cpd00242_e','EX_cpd00165_e','EX_cpd00418x_e'};
    community.Xfeed= [0 0];                                 
    community.Lfeed = [0 0 0 0 0 0 0 0 0 0 0 0];            
    community.LNMfeed = [0 0];                                  
    community.Gfeed = [0 0 0 0.21 0 0 0.78 0.0004];         
    community.Fin = 0;                                      
    community.Fout = 0;                                     
    community.Vgin = 0;                                     
    community.Vgout = 0;                                  
    
%Initialize variable bucket    
%   bucket: cumulative fluxes from selected species and internal cell reactions 
%       bucket.reac                         identifiers for internal reaction fluxes being monitored (mmol/L)
%       bucket.liquidreactormassbalance     metabolite species in liquid phase (mmol/L)
%       bucket.NMliquidreactormassbalance   non-metabolite species in liquid phase (mmol/L)
%       bucket.gasreactormassbalance        species in gas phase (mmol/L)
%       bucket.masstransfer                 species being transferred between phases (mmol/L)
%       bucket.numbbucket                   biotic buckets (cumulative flux) for selected species (mmol/L)
%       bucket.numabucket                   abiotic buckets for selected species (mmol/L)
%       bucket.flux                         internal reaction flux buckets for selected reactions (mmol/L)

    bucket.reac = {'rxn14202a','rxn14202b','rxn14202c','rxn14202d','rxn14202e','rxn14202f','rxn14202g','rxn01806mod','rxn05795mod','rxnP460a','rxnP460b','rxn10113','rxn10125mod','rxn08975','rxn10122mod','rxn10043mod','rxn00058','rxn10042','rxn12750','rxn14173a','rxn14173b','rxn00567a','rxn00567b','rxn00567c','rxn09008','rxn12750a','rxn12750b','rxn14010','rxn10811','Cytochrome_Exchange1','Cytochrome_Exchange2','Cytochrome_Exchange3','rxn12750mod'};
    bucket.liquidreactormassbalance = [0 0 0 0 0 0 0 0 0 0 0 0];
    bucket.NMliquidreactormassbalance =[0 0];
    bucket.gasreactormassbalance = [0 0 0 0 0 0 0 0];
    bucket.masstransfer = [0 0 0 0 0 0 0 0];
    bucket.numbbucket = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    bucket.numabucket = [0 0 0 0 0 0 0];
    bucket.flux = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    
%Initial conditions for system variables    
    initV = 0.005;                                                  %L
    initVg = 0.155;                                                 %L
    initX = [0.00016 0.0];                                   %gDCW
    initL = [0.0 0 0 0.001172 0.025 0 0 0 0.0 0.00005 0 0];         %mmol
    initLNM = [0 0];                                                %mmol
    initG = [0 0 0 1.308 0 0 4.86 0.0025];                          %mmol
    initAB=[0 0 0 0 0 0 0];                                         %mmol  
    initBB=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];       %mmol
    initR=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %mmol        

initialConditions = [initV initVg initX initL initLNM initG initAB initBB initR];

%Set simulation time and call solver
%Segment 1

%New Model Switch
    community.species{1,2} = changeRxnBounds(community.species{1,2},'rxn00567b',0,'u');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'rxn00567b',0,'l');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'rxn00567c',0,'u');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'rxn00567c',0,'l');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'Maintenance',8,'u');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'Maintenance',8,'l');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'Transporter',10000,'u');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'EX_cpd11416_e',0.016,'u');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'EX_cpd11416_e',0.0,'u');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxn05795mod',1.3,'u');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxn01806mod',0.4,'u');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxnP460a',0.05,'u');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxnP460b',0.05,'u');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxn14202g',0,'l');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxn14202g',0.04788,'u');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxn14202b',0,'l');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxn14202b',0,'u');    

    tspan = [0:0.01:0.25];
  
    [t, V, Vg, X, L, LNM, G, AB, BB, R] = Co4(community, bucket, tspan, initialConditions);    

ttot=[t];
  
Vtot=(V);
Vgtot=(Vg);
Xtot=(X);
Ltot=(L);
LNMtot=(LNM);
Gtot=(G);
ABtot=(AB);
BBtot=(BB);
Rtot=(R);

lastV = size(V);
lastVg = size(Vg);
lastX = size(X);
lastL = size(L);
lastLNM = size(LNM);
lastG = size(G);
lastAB = size(AB);
lastBB=size(BB);
lastR = size(R);

initV4 = V(lastV(1),1); %L
initVg4 = Vg(lastVg(1),1); %L
initX4 = [X(lastX(1),1) X(lastX(1),2)]; %gDCW
initL4 = [L(lastL(1),1) L(lastL(1),2) L(lastL(1),3) L(lastL(1),4) L(lastL(1),5) L(lastL(1),6) L(lastL(1),7) L(lastL(1),8) L(lastL(1),9) L(lastL(1),10) L(lastL(1),11) L(lastL(1),12)];   %mmol
initLNM4 = [LNM(lastLNM(1),1) LNM(lastLNM(1),2)]; %mmol
initG4 = [G(lastG(1),1) G(lastG(1),2) G(lastG(1),3) G(lastG(1),4) G(lastG(1),5) G(lastG(1),6) G(lastG(1),7) G(lastG(1),8)];%mmol
initAB4 = [AB(lastAB(1),1) AB(lastAB(1),2) AB(lastAB(1),3) AB(lastAB(1),4) AB(lastAB(1),5) AB(lastAB(1),6) AB(lastAB(1),7)]; %mmol/L
initBB4 = [BB(lastBB(1),1) BB(lastBB(1),2) BB(lastBB(1),3) BB(lastBB(1),4) BB(lastBB(1),5) BB(lastBB(1),6) BB(lastBB(1),7) BB(lastBB(1),8) BB(lastBB(1),9) BB(lastBB(1),10) BB(lastBB(1),11) BB(lastBB(1),12) BB(lastBB(1),13) BB(lastBB(1),14) BB(lastBB(1),15) BB(lastBB(1),16) BB(lastBB(1),17) BB(lastBB(1),18) BB(lastBB(1),19) BB(lastBB(1),20) BB(lastBB(1),21) BB(lastBB(1),22)  BB(lastBB(1),23) BB(lastBB(1),24)]; %mmol/L
initR4 = [R(lastR(1),1) R(lastR(1),2) R(lastR(1),3) R(lastR(1),4) R(lastR(1),5) R(lastR(1),6) R(lastR(1),7) R(lastR(1),8) R(lastR(1),9) R(lastR(1),10) R(lastR(1),11) R(lastR(1),12) R(lastR(1),13) R(lastR(1),14) R(lastR(1),15) R(lastR(1),16) R(lastR(1),17) R(lastR(1),18) R(lastR(1),19) R(lastR(1),20) R(lastR(1),21) R(lastR(1),22) R(lastR(1),23) R(lastR(1),24) R(lastR(1),25) R(lastR(1),26) R(lastR(1),27) R(lastR(1),28) R(lastR(1),29) R(lastR(1),30) R(lastR(1),31) R(lastR(1),32) R(lastR(1),33) R(lastR(1),34) R(lastR(1),35) R(lastR(1),36) R(lastR(1),37) R(lastR(1),38) R(lastR(1),39) R(lastR(1),40) R(lastR(1),41) R(lastR(1),42) R(lastR(1),43) R(lastR(1),44) R(lastR(1),45) R(lastR(1),46) R(lastR(1),47) R(lastR(1),48) R(lastR(1),49) R(lastR(1),50) R(lastR(1),51) R(lastR(1),52) R(lastR(1),53) R(lastR(1),54) R(lastR(1),55) R(lastR(1),56) R(lastR(1),57) R(lastR(1),58) R(lastR(1),59) R(lastR(1),60) R(lastR(1),61) R(lastR(1),62) R(lastR(1),63) R(lastR(1),64) R(lastR(1),65) R(lastR(1),66)]; %mmol/L

initialConditions = [initV4 initVg4 initX4 initL4 initLNM4 initG4 initAB4 initBB4 initR4];

%Set simulation time and call solver
%Segment 1A

%New Model Switch
    community.species{1,2} = changeRxnBounds(community.species{1,2},'rxn00567b',0,'u');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'rxn00567b',0,'l');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'rxn00567c',0,'u');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'rxn00567c',0,'l');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'Maintenance',8,'u');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'Maintenance',8,'l');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'Transporter',10000,'u');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'EX_cpd11416_e',0.016,'u');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'EX_cpd11416_e',0.0,'u');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxn05795mod',1.3,'u');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxn01806mod',0.4,'u');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxnP460a',0.05,'u');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxnP460b',0.05,'u');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxn14202g',0,'l');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxn14202g',0.04788,'u');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxn14202b',0,'l');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxn14202b',0,'u');    

    tspan = [0.25:0.01:2];
  
    [t, V, Vg, X, L, LNM, G, AB, BB, R] = Co(community, bucket, tspan, initialConditions);    

ttot=[t];
  
ttot=[ttot; t];
Vtot=[Vtot; V];
Vgtot=[Vgtot; Vg];
Xtot=[Xtot; X];
Ltot=[Ltot; L];
LNMtot = [LNMtot;LNM];
Gtot=[Gtot; G];
ABtot=[ABtot;AB];
BBtot=[BBtot; BB];
Rtot=[Rtot; R];

lastV = size(V);
lastVg = size(Vg);
lastX = size(X);
lastL = size(L);
lastLNM = size(LNM);
lastG = size(G);
lastAB = size(AB);
lastBB=size(BB);
lastR = size(R);

initV1 = V(lastV(1),1); %L
initVg1 = Vg(lastVg(1),1); %L
initX1 = [X(lastX(1),1) X(lastX(1),2)]; %gDCW
initL1 = [L(lastL(1),1) L(lastL(1),2) L(lastL(1),3) L(lastL(1),4) L(lastL(1),5) L(lastL(1),6) L(lastL(1),7) L(lastL(1),8) L(lastL(1),9) L(lastL(1),10) L(lastL(1),11) L(lastL(1),12)];   %mmol
initLNM1 = [LNM(lastLNM(1),1) LNM(lastLNM(1),2)]; %mmol
initG1 = [G(lastG(1),1) G(lastG(1),2) G(lastG(1),3) G(lastG(1),4) G(lastG(1),5) G(lastG(1),6) G(lastG(1),7) G(lastG(1),8)];%mmol
initAB1 = [AB(lastAB(1),1) AB(lastAB(1),2) AB(lastAB(1),3) AB(lastAB(1),4) AB(lastAB(1),5) AB(lastAB(1),6) AB(lastAB(1),7)]; %mmol/L
initBB1 = [BB(lastBB(1),1) BB(lastBB(1),2) BB(lastBB(1),3) BB(lastBB(1),4) BB(lastBB(1),5) BB(lastBB(1),6) BB(lastBB(1),7) BB(lastBB(1),8) BB(lastBB(1),9) BB(lastBB(1),10) BB(lastBB(1),11) BB(lastBB(1),12) BB(lastBB(1),13) BB(lastBB(1),14) BB(lastBB(1),15) BB(lastBB(1),16) BB(lastBB(1),17) BB(lastBB(1),18) BB(lastBB(1),19) BB(lastBB(1),20) BB(lastBB(1),21) BB(lastBB(1),22)  BB(lastBB(1),23) BB(lastBB(1),24)]; %mmol/L
initR1 = [R(lastR(1),1) R(lastR(1),2) R(lastR(1),3) R(lastR(1),4) R(lastR(1),5) R(lastR(1),6) R(lastR(1),7) R(lastR(1),8) R(lastR(1),9) R(lastR(1),10) R(lastR(1),11) R(lastR(1),12) R(lastR(1),13) R(lastR(1),14) R(lastR(1),15) R(lastR(1),16) R(lastR(1),17) R(lastR(1),18) R(lastR(1),19) R(lastR(1),20) R(lastR(1),21) R(lastR(1),22) R(lastR(1),23) R(lastR(1),24) R(lastR(1),25) R(lastR(1),26) R(lastR(1),27) R(lastR(1),28) R(lastR(1),29) R(lastR(1),30) R(lastR(1),31) R(lastR(1),32) R(lastR(1),33) R(lastR(1),34) R(lastR(1),35) R(lastR(1),36) R(lastR(1),37) R(lastR(1),38) R(lastR(1),39) R(lastR(1),40) R(lastR(1),41) R(lastR(1),42) R(lastR(1),43) R(lastR(1),44) R(lastR(1),45) R(lastR(1),46) R(lastR(1),47) R(lastR(1),48) R(lastR(1),49) R(lastR(1),50) R(lastR(1),51) R(lastR(1),52) R(lastR(1),53) R(lastR(1),54) R(lastR(1),55) R(lastR(1),56) R(lastR(1),57) R(lastR(1),58) R(lastR(1),59) R(lastR(1),60) R(lastR(1),61) R(lastR(1),62) R(lastR(1),63) R(lastR(1),64) R(lastR(1),65) R(lastR(1),66)]; %mmol/L

initialConditions = [initV1 initVg1 initX1 initL1 initLNM1 initG1 initAB1 initBB1 initR1];

%Segment 2

%Reduction in Nitrous and Nitric Oxide Production
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxn05795mod',1.3,'u');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxn01806mod',0.4,'u');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxnP460a',0.05,'u');
    community.species{1,1} = changeRxnBounds(community.species{1,1},'rxnP460b',0.05,'u');

tspan = [2:0.01:3];       

    [t, V, Vg, X, L, LNM, G, AB, BB, R] = Co2(community, bucket, tspan, initialConditions);    

ttot=[t];
  
ttot=[ttot; t];
Vtot=[Vtot; V];
Vgtot=[Vgtot; Vg];
Xtot=[Xtot; X];
Ltot=[Ltot; L];
LNMtot = [LNMtot;LNM];
Gtot=[Gtot; G];
ABtot=[ABtot;AB];
BBtot=[BBtot; BB];
Rtot=[Rtot; R];

lastV = size(V);
lastVg = size(Vg);
lastX = size(X);
lastL = size(L);
lastLNM = size(LNM);
lastG = size(G);
lastAB = size(AB);
lastBB=size(BB);
lastR = size(R);

initV2 = V(lastV(1),1); %L
initVg2 = Vg(lastVg(1),1); %L
initX2 = [X(lastX(1),1) X(lastX(1),2)]; %gDCW
initL2 = [L(lastL(1),1) L(lastL(1),2) L(lastL(1),3) L(lastL(1),4) L(lastL(1),5) L(lastL(1),6) L(lastL(1),7) L(lastL(1),8) L(lastL(1),9) L(lastL(1),10) L(lastL(1),11) L(lastL(1),12)];   %mmol
initLNM2 = [LNM(lastLNM(1),1) LNM(lastLNM(1),2)]; %mmol
initG2= [G(lastG(1),1) G(lastG(1),2) G(lastG(1),3) G(lastG(1),4) G(lastG(1),5) G(lastG(1),6) G(lastG(1),7) G(lastG(1),8)];%mmol
initAB2 = [AB(lastAB(1),1) AB(lastAB(1),2) AB(lastAB(1),3) AB(lastAB(1),4) AB(lastAB(1),5) AB(lastAB(1),6) AB(lastAB(1),7)]; %mmol/L
initBB2 = [BB(lastBB(1),1) BB(lastBB(1),2) BB(lastBB(1),3) BB(lastBB(1),4) BB(lastBB(1),5) BB(lastBB(1),6) BB(lastBB(1),7) BB(lastBB(1),8) BB(lastBB(1),9) BB(lastBB(1),10) BB(lastBB(1),11) BB(lastBB(1),12) BB(lastBB(1),13) BB(lastBB(1),14) BB(lastBB(1),15) BB(lastBB(1),16) BB(lastBB(1),17) BB(lastBB(1),18) BB(lastBB(1),19) BB(lastBB(1),20) BB(lastBB(1),21) BB(lastBB(1),22)  BB(lastBB(1),23) BB(lastBB(1),24)]; %mmol/L
initR2 = [R(lastR(1),1) R(lastR(1),2) R(lastR(1),3) R(lastR(1),4) R(lastR(1),5) R(lastR(1),6) R(lastR(1),7) R(lastR(1),8) R(lastR(1),9) R(lastR(1),10) R(lastR(1),11) R(lastR(1),12) R(lastR(1),13) R(lastR(1),14) R(lastR(1),15) R(lastR(1),16) R(lastR(1),17) R(lastR(1),18) R(lastR(1),19) R(lastR(1),20) R(lastR(1),21) R(lastR(1),22) R(lastR(1),23) R(lastR(1),24) R(lastR(1),25) R(lastR(1),26) R(lastR(1),27) R(lastR(1),28) R(lastR(1),29) R(lastR(1),30) R(lastR(1),31) R(lastR(1),32) R(lastR(1),33) R(lastR(1),34) R(lastR(1),35) R(lastR(1),36) R(lastR(1),37) R(lastR(1),38) R(lastR(1),39) R(lastR(1),40) R(lastR(1),41) R(lastR(1),42) R(lastR(1),43) R(lastR(1),44) R(lastR(1),45) R(lastR(1),46) R(lastR(1),47) R(lastR(1),48) R(lastR(1),49) R(lastR(1),50) R(lastR(1),51) R(lastR(1),52) R(lastR(1),53) R(lastR(1),54) R(lastR(1),55) R(lastR(1),56) R(lastR(1),57) R(lastR(1),58) R(lastR(1),59) R(lastR(1),60) R(lastR(1),61) R(lastR(1),62) R(lastR(1),63) R(lastR(1),64) R(lastR(1),65) R(lastR(1),66)]; %mmol/L

initialConditions = [initV2 initVg2 initX2 initL2 initLNM2 initG2 initAB2 initBB2 initR2];

%Segment 3

%Poughon Model Switch
    community.species{1,2} = changeRxnBounds(community.species{1,2},'rxn00567b',1.3,'u');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'rxn00567b',0,'l');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'rxn00567c',1.3,'u');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'rxn00567c',-1.3,'l');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'Maintenance',18.52,'u');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'Maintenance',18.52,'l');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'Transporter',0,'u');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'rxn12750b',10000,'u');
    community.species{1,2} = changeRxnBounds(community.species{1,2},'rxn12750b',-10000,'l');
    
tspan = [3:0.01:4];
    
    [t, V, Vg, X, L, LNM, G, AB, BB, R] = Co3(community, bucket, tspan, initialConditions);       

ttot=[ttot; t];
Vtot=[Vtot; V];
Vgtot=[Vgtot; Vg];
Xtot=[Xtot; X];
Ltot=[Ltot; L];
LNMtot = [LNMtot;LNM];
Gtot=[Gtot; G];
ABtot=[ABtot;AB];
BBtot=[BBtot; BB];
Rtot=[Rtot; R];

lastV = size(V);
lastVg = size(Vg);
lastX = size(X);
lastL = size(L);
lastLNM = size(LNM);
lastG = size(G);
lastAB = size(AB);
lastBB=size(BB);
lastR = size(R);

initV3 = V(lastV(1),1); %L
initVg3 = Vg(lastVg(1),1); %L
initX3 = [X(lastX(1),1) X(lastX(1),2)]; %gDCW
initL3 = [L(lastL(1),1) L(lastL(1),2) L(lastL(1),3) L(lastL(1),4) L(lastL(1),5) L(lastL(1),6) L(lastL(1),7) L(lastL(1),8) L(lastL(1),9) L(lastL(1),10) L(lastL(1),11) L(lastL(1),12)];   %mmol
initLNM3 = [LNM(lastLNM(1),1) LNM(lastLNM(1),2)]; %mmol
initG3= [G(lastG(1),1) G(lastG(1),2) G(lastG(1),3) G(lastG(1),4) G(lastG(1),5) G(lastG(1),6) G(lastG(1),7) G(lastG(1),8)];%mmol
initAB3 = [AB(lastAB(1),1) AB(lastAB(1),2) AB(lastAB(1),3) AB(lastAB(1),4) AB(lastAB(1),5) AB(lastAB(1),6) AB(lastAB(1),7)]; %mmol/L
initBB3 = [BB(lastBB(1),1) BB(lastBB(1),2) BB(lastBB(1),3) BB(lastBB(1),4) BB(lastBB(1),5) BB(lastBB(1),6) BB(lastBB(1),7) BB(lastBB(1),8) BB(lastBB(1),9) BB(lastBB(1),10) BB(lastBB(1),11) BB(lastBB(1),12) BB(lastBB(1),13) BB(lastBB(1),14) BB(lastBB(1),15) BB(lastBB(1),16) BB(lastBB(1),17) BB(lastBB(1),18) BB(lastBB(1),19) BB(lastBB(1),20) BB(lastBB(1),21) BB(lastBB(1),22)  BB(lastBB(1),23) BB(lastBB(1),24)]; %mmol/L
initR3 = [R(lastR(1),1) R(lastR(1),2) R(lastR(1),3) R(lastR(1),4) R(lastR(1),5) R(lastR(1),6) R(lastR(1),7) R(lastR(1),8) R(lastR(1),9) R(lastR(1),10) R(lastR(1),11) R(lastR(1),12) R(lastR(1),13) R(lastR(1),14) R(lastR(1),15) R(lastR(1),16) R(lastR(1),17) R(lastR(1),18) R(lastR(1),19) R(lastR(1),20) R(lastR(1),21) R(lastR(1),22) R(lastR(1),23) R(lastR(1),24) R(lastR(1),25) R(lastR(1),26) R(lastR(1),27) R(lastR(1),28) R(lastR(1),29) R(lastR(1),30) R(lastR(1),31) R(lastR(1),32) R(lastR(1),33) R(lastR(1),34) R(lastR(1),35) R(lastR(1),36) R(lastR(1),37) R(lastR(1),38) R(lastR(1),39) R(lastR(1),40) R(lastR(1),41) R(lastR(1),42) R(lastR(1),43) R(lastR(1),44) R(lastR(1),45) R(lastR(1),46) R(lastR(1),47) R(lastR(1),48) R(lastR(1),49) R(lastR(1),50) R(lastR(1),51) R(lastR(1),52) R(lastR(1),53) R(lastR(1),54) R(lastR(1),55) R(lastR(1),56) R(lastR(1),57) R(lastR(1),58) R(lastR(1),59) R(lastR(1),60) R(lastR(1),61) R(lastR(1),62) R(lastR(1),63) R(lastR(1),64) R(lastR(1),65) R(lastR(1),66)]; %mmol/L

initialConditions = [initV3 initVg3 initX3 initL3 initLNM3 initG3 initAB3 initBB3 initR3];

%Segment 4

tspan = [4:0.01:10];

    [t, V, Vg, X, L, LNM, G, AB, BB, R] = Co1(community, bucket, tspan, initialConditions);       

ttot=[ttot; t];
Vtot=[Vtot; V];
Vgtot=[Vgtot; Vg];
Xtot=[Xtot; X];
Ltot=[Ltot; L];
LNMtot = [LNMtot;LNM];
Gtot=[Gtot; G];
ABtot=[ABtot;AB];
BBtot=[BBtot; BB];
Rtot=[Rtot; R];
