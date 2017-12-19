%Rate constants and other parameters
ph=-7.4;kla1 = 32.67; kla2 = 35.85; kla3 = 0.000956; kla4 = 1.982; kla5 = 70.71; kla6 = 1.391; kga1 = 19.764; T=303.15; kga2=7416; KsNO = 1.3*10^-2;KsNO2 = 1.3*10^-2;
k1 = 10.8*10^10 ; kminus1=1.8*10^11; k2 = 8640; k3 = 24.84; k4 = 102.96; kminus4 = 6.012*10^8; k5 = 279360000;
k6=133.2;kminus6=96.12;k7=7.6*10^-6;

%Liquid Volume (L)
Vdum(1)=y(1);

%Gas volume (L)
Vdum(2)=y(2);

%Biomass Concentration Transformation (gDCW/L)
Xdum=y(3:2+numSpecies)/y(1);

%Liquid Phase Species Concentration Transformation (mmol/L)
Ldum=y(2+numSpecies+1:2+numSpecies+numLiquidSpecies)/y(1);

%Liquid Phase Non-Metabolite Species Concentration Transformation (mmol/L)
LNMdum=y(2+numSpecies+numLiquidSpecies+1:2+numSpecies+numLiquidSpecies+numLiquidNMSpecies)/y(1);

%Gas Phase Species Concentration Transformation (atm)
Gdum=y(2+numSpecies+numLiquidSpecies+numLiquidNMSpecies+1:2+numSpecies+numLiquidSpecies+numLiquidNMSpecies+numGasSpecies)*(0.00008206*T)/y(2);

%Liquid species reactor balance (mmol)
LBdum=((community.Fin/y(1))*(community.Lfeed' - Ldum)*y(1));

%Liquid Non-Metabolite Species reactor balance (mmol)
LNMBdum=(community.Fin/y(1)*(community.LNMfeed' - LNMdum)*y(1));

%Gas species reactor balance (mmol)
GBdum=(community.Vgin/y(2)*(community.Gfeed' - Gdum)*(y(2)/0.00008206*T));

%Biotic Reaction Buckets (mmol/(L))
BBdum(1)= Xdum(1)*vs(1,1); %NO2- B(1)
BBdum(2)= Xdum(2)*vs(2,1); %NO2- B(2)
BBdum(3)= Xdum(1)*vs(1,2); %nitrate B(3)
BBdum(4)= Xdum(2)*vs(2,2); %nitrate B(4)
BBdum(5)= Xdum(1)*vs(1,3); %NO B(5)
BBdum(6)= Xdum(2)*vs(2,3); %NO B(6)
BBdum(7)= Xdum(1)*vs(1,4); %O2 B(7)
BBdum(8)= Xdum(2)*vs(2,4); %O2 B(8)
BBdum(9)= Xdum(1)*vs(1,5); %NH3 B(9)
BBdum(10)= Xdum(2)*vs(2,5); %NH3 B(10)
BBdum(11)= Xdum(1)*vs(1,6); %CHO2 B(11)
BBdum(12)= Xdum(2)*vs(2,6); %CHO2 B(12)
BBdum(13)= Xdum(1)*vs(1,7); %N2O B(13)
BBdum(14)= Xdum(2)*vs(2,7); %N2O B(14)
BBdum(15)= Xdum(1)*vs(1,8); %N2 B(15)
BBdum(16)= Xdum(2)*vs(2,8); %N2 B(16)
BBdum(17)= Xdum(1)*vs(1,9); %CO2 B(17)
BBdum(18)= Xdum(2)*vs(2,9); %CO2 B(18)
BBdum(19)= Xdum(1)*vs(1,10); %HCO3 B(19)
BBdum(20)= Xdum(2)*vs(2,10); %HCO3 B(20)
BBdum(21)= Xdum(1)*vs(1,11); %NH2OH B(21)
BBdum(22)= Xdum(2)*vs(2,11); %NH2OH B(22)
BBdum(23)= Xdum(1)*vs(1,12); %NOhao B(23)
BBdum(24)= Xdum(2)*vs(2,12); %NOhao B(24)

%Mass transfer Buckets(mmol/(L))
MBdum(1)= kla1*(1.269*Gdum(1)-Ldum(3)); %NO
MBdum(2)= -kga1*(Gdum(2)-0.041667*LNMdum(1))*y(2)/(0.00008206*T); %NO2
MBdum(3)= -kga2*(Gdum(3)-0.00002*LNMdum(2))*y(2)/(0.00008206*T); %HONO
MBdum(4)= kla2*(1.15600*Gdum(4)-Ldum(4)); %O2
MBdum(5)= kla3*(43374*Gdum(5)-Ldum(5)); %NH3
MBdum(6)= kla4*(20.91*Gdum(6)-Ldum(7)); %N2O
MBdum(7)= kla5*(0.586*Gdum(7)-Ldum(8)); %N2
MBdum(8)= kla6*(29.78*Gdum(8)-Ldum(9)); %CO2

%Abiotic reaction buckets (mmol/(L))
ABdum(1)= (k1*LNMdum(2)-kminus1*Ldum(1)*1000*(10^(ph))); %liquid formation of NO2- from HONO A(1)
ABdum(2)= 4*k2*((Ldum(3)))*((Ldum(3)))*Ldum(4); %liquid oxidation of 4 NO to 4 NO2- A(2)
ABdum(3)= 2*k3*Gdum(1)*Gdum(1)*Gdum(4)/((0.00008206*T)^3); %gas oxidation of 2 NO2 to 2 NO A(3)
ABdum(4)= (k4*LNMdum(2)*LNMdum(2)-kminus4*((Ldum(3)))*LNMdum(1)); %liquid side oxidation of 2 HONO to NO and NO2 A(4)
ABdum(5)= (k5*LNMdum(1)*LNMdum(1)); %liquid side oxidation of 2NO2 to NO3- and NO2- A(5)
ABdum(6)= -k6*Ldum(9)+kminus6*1000*(10^ph)*Ldum(10); %liquid side carbonate equilibrium A(6)
ABdum(7)= k7*Ldum(11)*1000*(10^ph)*LNMdum(2); %liquid side formation of N2O from nitrite and hydroxlamine A(7)

%Bioreactor volume changes (L)
dy(1) = community.Fin - community.Fout;
dy(2) = community.Vgin - community.Vgout;

%Biomass changes (gDCW)
dy(3)= y(1)*(mu(1)*Xdum(1)+community.Fin/y(1)*(community.Xfeed(1) - Xdum(1))); %Neuro Biomass X(1)
dy(4)= y(1)*(mu(2)*Xdum(2)+community.Fin/y(1)*(community.Xfeed(2) - Xdum(2))); %Nwino Biomass X(2)
   
%Liquid phase metabolite reaction network (mmol)
dy(5)= y(1)*(ABdum(1)+ABdum(2)+BBdum(1)+BBdum(2)-ABdum(7))+LBdum(1); %NO2- L(1)
dy(6)= y(1)*(ABdum(5)+BBdum(3)+BBdum(4))+LBdum(2); %nitrate L(2)
dy(7)= y(1)*MBdum(1)+y(1)*(ABdum(4)- ABdum(2)+BBdum(5)+BBdum(23)+BBdum(6)+BBdum(24))+LBdum(3); %NO L(3)
dy(8)= y(1)*MBdum(4)+y(1)*(-0.25*ABdum(2)+BBdum(7)+BBdum(8))+LBdum(4); %O2 L(4)
dy(9)= y(1)*MBdum(5)+y(1)*(BBdum(9)+ BBdum(10))+LBdum(5); %NH3 L(5)
dy(10)= y(1)*(BBdum(11)+ BBdum(12))+LBdum(6); %CHO2 L(6)
dy(11)= y(1)*MBdum(6)+y(1)*(BBdum(13)+ BBdum(14)+ABdum(7))+LBdum(7); %N2O L(7)
dy(12)= y(1)*MBdum(7)+y(1)*(BBdum(15)+ BBdum(16))+LBdum(8); %N2 L(8)
dy(13)= y(1)*MBdum(8)+y(1)*(ABdum(6)+BBdum(17)+BBdum(18))+LBdum(9); %CO2 L(9)
dy(14)= y(1)*(-ABdum(6)+BBdum(19)+ BBdum(20))+LBdum(10); %HCO3 L(10)
dy(15)=y(1)*(BBdum(21)+BBdum(22))+LBdum(11); %L(11) NH2OH
dy(16)=y(1)*(BBdum(23)+BBdum(24))+LBdum(12); %L(12) NOHAO

%Liquid phase non-metabolite reaction network (mmol)
dy(17)= y(1)*(ABdum(4)-2*ABdum(5))-y(2)*MBdum(2)+LNMBdum(1); %NO2 LNM(1)
dy(18)= y(1)*(-ABdum(1)-2*ABdum(4)+ABdum(5))-y(2)*MBdum(3)+LNMBdum(2);%HONO LNM(2)

%Gas phase reaction/mass transfer network (mmol)
dy(19)= -y(1)*MBdum(1)-y(2)*ABdum(3)+GBdum(1)*y(2)/(0.00008206*T); %NO(g)G(1)
dy(20)= y(2)*MBdum(2)+y(2)*ABdum(3)+GBdum(2)*y(2)/(0.00008206*T);%NO2(g) G(2)
dy(21)= y(2)*MBdum(3)+ GBdum(3)*y(2)/(0.00008206*T); %HONO (g) G(3)
dy(22)= -y(1)*MBdum(4)-0.5*y(2)*ABdum(3)+GBdum(4)*y(2)/(0.00008206*T); %O2 (g) G(4)
dy(23)= -y(1)*MBdum(5)+GBdum(5)*y(2)/(0.00008206*T); %NH3 (g) G(5)
dy(24)= -y(1)*MBdum(6)+GBdum(6)*y(2)/(0.00008206*T); %N2O (g) G(6)
dy(25)= -y(1)*MBdum(7)+GBdum(7)*y(2)/(0.00008206*T); %N2 (g) G(7)
dy(26)= -y(1)*MBdum(8)+GBdum(8)*y(2)/(0.00008206*T); %CO2 (g) G(8)

%Abiotic Buckets Retrieval (mmol/(L))
dy(27)= ABdum(1);
dy(28)= ABdum(2);
dy(29)= ABdum(3);
dy(30)= ABdum(4);
dy(31)= ABdum(5);
dy(32)= ABdum(6);
dy(33)= ABdum(7);

%Biotic Buckets Retrieval (mmol/(L))
dy(34)= Xdum(1)*vs(1,1); %NO2- B(1)
dy(35)= Xdum(2)*vs(2,1); %NO2- B(2)
dy(36)= Xdum(1)*vs(1,2); %nitrate B(3)
dy(37)= Xdum(2)*vs(2,2); %nitrate B(4)
dy(38)= Xdum(1)*vs(1,3); %NO B(5)
dy(39)= Xdum(2)*vs(2,3); %NO B(6)
dy(40)= Xdum(1)*vs(1,4); %O2 B(7)
dy(41)= Xdum(2)*vs(2,4); %O2 B(8)
dy(42)= Xdum(1)*vs(1,5); %NH3 B(9)
dy(43)= Xdum(2)*vs(2,5); %NH3 B(10)
dy(44)= Xdum(1)*vs(1,6); %CHO2 B(11)
dy(45)= Xdum(2)*vs(2,6); %CHO2 B(12)
dy(46)= Xdum(1)*vs(1,7); %N2O B(13)
dy(47)= Xdum(2)*vs(2,7); %N2O B(14)
dy(48)= Xdum(1)*vs(1,8); %N2 B(15)
dy(49)= Xdum(2)*vs(2,8); %N2 B(16)
dy(50)= Xdum(1)*vs(1,9); %CO2 B(17)
dy(51)= Xdum(2)*vs(2,9); %CO2 B(18)
dy(52)= Xdum(1)*vs(1,10); %HCO3 B(19)
dy(53)= Xdum(2)*vs(2,10); %HCO3 B(20)
dy(54)= Xdum(1)*vs(1,11); %NH2OH B(21)
dy(55)= Xdum(2)*vs(2,11); %NH2OH B(22)
dy(56)= Xdum(1)*vs(1,12); %NOhao B(23)
dy(57)= Xdum(2)*vs(2,12); %NOhao B(24)

%Internal Flux Retrieval (mmol/(L))
dy(58)=Xdum(1)*vs(1,13); %R(1) rxn14202a
dy(59)=Xdum(2)*vs(2,13); %R(1) 
dy(60)=Xdum(1)*vs(1,14); %R(2) rxn14202b
dy(61)=Xdum(2)*vs(2,14); %R(2) 
dy(62)=Xdum(1)*vs(1,15); %R(3) rxn14202c
dy(63)=Xdum(2)*vs(2,15); %R(3) 
dy(64)=Xdum(1)*vs(1,16); %R(4) rxn14202d
dy(65)=Xdum(2)*vs(2,16); %R(4) 
dy(66)=Xdum(1)*vs(1,17); %R(5) rxn014202e
dy(67)=Xdum(2)*vs(2,17); %R(5) 
dy(68)=Xdum(1)*vs(1,18); %R(6) rxn14202f
dy(69)=Xdum(2)*vs(2,18); %R(6)
dy(70)=Xdum(1)*vs(1,19); %R(7) rxn14202g
dy(71)=Xdum(2)*vs(2,19); %R(7) 
dy(72)=Xdum(1)*vs(1,20); %R(8) rxn1806mod
dy(73)=Xdum(2)*vs(2,20); %R(8) 
dy(74)=Xdum(1)*vs(1,21); %R(9) rxn05795mod
dy(75)=Xdum(2)*vs(2,21); %R(9) 
dy(76)=Xdum(1)*vs(1,22); %R(10) rxnP460a
dy(77)=Xdum(2)*vs(2,22); %R(10) 
dy(78)=Xdum(1)*vs(1,23); %R(11) rxnP460b
dy(79)=Xdum(2)*vs(2,23); %R(11) 
dy(80)=Xdum(1)*vs(1,24); %R(12) rxn10113
dy(81)=Xdum(2)*vs(2,24); %R(12) 
dy(82)=Xdum(1)*vs(1,25); %R(13) rxn10125mod
dy(83)=Xdum(2)*vs(2,25); %R(13)
dy(84)=Xdum(1)*vs(1,26); %R(14) rxn08975
dy(85)=Xdum(2)*vs(2,26); %R(14) 
dy(86)=Xdum(1)*vs(1,27); %R(15) rxn10122mod
dy(87)=Xdum(2)*vs(2,27); %R(15) 
dy(88)=Xdum(1)*vs(1,28); %R(16) rxn10043mod
dy(89)=Xdum(2)*vs(2,28); %R(16) 
dy(90)=Xdum(1)*vs(1,29); %R(17) rxn00058
dy(91)=Xdum(2)*vs(2,29); %R(17) 
dy(92)=Xdum(1)*vs(1,30); %R(18) rxn10042
dy(93)=Xdum(2)*vs(2,30); %R(18) 
dy(94)=Xdum(1)*vs(1,31); %R(19) rxn12750a
dy(95)=Xdum(2)*vs(2,31); %R(19) 
dy(96)=Xdum(1)*vs(1,32); %R(20) rxn14173a
dy(97)=Xdum(2)*vs(2,32); %R(20)
dy(98)=Xdum(1)*vs(1,33); %R(21) rxn14173b
dy(99)=Xdum(2)*vs(2,33); %R(21) 
dy(100)=Xdum(1)*vs(1,34); %R(22) rxn00567a
dy(101)=Xdum(2)*vs(2,34); %R(22) 
dy(102)=Xdum(1)*vs(1,35); %R(23) rxn00567b
dy(103)=Xdum(2)*vs(2,35); %R(23) 
dy(104)=Xdum(1)*vs(1,36); %R(24) rxn00567c
dy(105)=Xdum(2)*vs(2,36); %R(24) 
dy(106)=Xdum(1)*vs(1,37); %R(25) rxn09008
dy(107)=Xdum(2)*vs(2,37); %R(25) 
dy(108)=Xdum(1)*vs(1,38); %R(26) rxn12750a
dy(109)=Xdum(2)*vs(2,38); %R(26) 
dy(110)=Xdum(1)*vs(1,39); %R(27) rxn12750b
dy(111)=Xdum(2)*vs(2,39); %R(27) 
dy(112)=Xdum(1)*vs(1,40); %R(28) rxn14010
dy(113)=Xdum(2)*vs(2,40); %R(28) 
dy(114)=Xdum(1)*vs(1,41); %R(29) rxn10811
dy(115)=Xdum(2)*vs(2,41); %R(29) 
dy(116)=Xdum(1)*vs(1,42); %R(30) Cytochrome Exchange1
dy(117)=Xdum(2)*vs(2,42); %R(30) 
dy(118)=Xdum(1)*vs(1,43); %R(31) Cytochrome Exchange2
dy(119)=Xdum(2)*vs(2,43); %R(31) 
dy(120)=Xdum(1)*vs(1,44); %R(32) Cytochrome Exchange3
dy(121)=Xdum(2)*vs(2,44); %R(32) 
dy(122)=Xdum(1)*vs(1,45); %R(33) 12750b
dy(123)=Xdum(2)*vs(2,45); %R(33) 