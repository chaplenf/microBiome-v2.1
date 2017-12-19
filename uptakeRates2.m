% Updates the uptake rates of the organisms based on liquid concentrations
% This code is called in coCulture_ode prior to FBA calculation
for i = 1:numSpecies
    switch community.species{1,i}.description
        case 'Nitrosomonas europaea'			
            % ammonia
            [~, loc] = ismember('EX_cpd00013_e',community.mets);
            C = subFlux(loc)/V;
            Vmax_amml = -61.34*C/(C+0.0365); 
            Vmax_ammu = 10000;
            community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00013_e',[Vmax_amml],'l');
			community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00013_e',[Vmax_ammu],'u');			
			% nitrite
            [~, loc] = ismember('EX_cpd00075_e',community.mets);
            C = subFlux(loc)/V;
            Vmax_NO2l = -61.34*C/(C+0.0365); 
            Vmax_NO2u=  10000;
            community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00075_e',[Vmax_NO2l],'l');
            community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00075_e',[Vmax_NO2u],'u');      
            % CO2
            [TF, loc] = ismember('EX_cpd00011_e',community.mets);
            C = subFlux(loc)/V;
            Vmax_CO2l = -61.34*C/(C+0.0038);  
            Vmax_CO2u = 10000;
     		community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00011_e',[Vmax_CO2l],'l');
			community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00011_e',[Vmax_CO2u],'u');			
			% bicarb
            [~, loc] = ismember('EX_cpd00242_e',community.mets);
            C = subFlux(loc)/V;
            Vmax_HCO3l = -61.34*C/(C+0.0038);
            Vmax_HCO3u=  10000;
            community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00242_e',[Vmax_HCO3l],'l');
            community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00242_e',[Vmax_HCO3u],'u');      
            % NO
            [~, loc] = ismember('EX_cpd00418_e',community.mets);
            C = subFlux(loc)/V;
            Vmax_NOl =-61.34*C/(C+0.01275);
            Vmax_NOu=10000;
			community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00418_e',[Vmax_NOl],'l');
			community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00418_e',[Vmax_NOu],'u');  
            % NOx
            [~, loc] = ismember('EX_cpd00418_e',community.mets);
            C = subFlux(loc)/V;
            Vmax_NOl=-61.34*C/(C+0.01275);
            Vmax_NOu=10000;
			community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00418x_c',[Vmax_NOl],'l');
			community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00418x_c',[Vmax_NOu],'u');  
            % NO3
            [~, loc] = ismember('EX_cpd00209_e',community.mets);
            C = subFlux(loc)/V;
            Vmax_NO3 = 0;
			community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00209_e',[Vmax_NO3],'l');     
            % O2
            [~, loc] = ismember('EX_cpd00007_e',community.mets);
            C = subFlux(loc)/V;
            Vmax_O2l = -92.01*C/(C+0.01275);
            Vmax_O2u = 10000;
			community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00007_e',[Vmax_O2l],'l');
            community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00007_e',[Vmax_O2u],'u');
            % CHO2
           [~, loc] = ismember('EX_cpd00047_e',community.mets);
            C = subFlux(loc)/V;
            Vmax_CHO2l = 0;  
            Vmax_CHO2u = 10000;
            community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00047_e',[Vmax_CHO2l],'l');
            community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00047_e',[Vmax_CHO2u],'u');
            % N2O
            [~, loc] = ismember('EX_cpd00659_e',community.mets);
            C = subFlux(loc)/V;
            Vmax_N2Ol = 0;  
            Vmax_N2Ou = 10000;
			community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00659_e',[Vmax_N2Ol],'l');
            community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00659_e',[Vmax_N2Ou],'u');                   
            % N2
            [~, loc] = ismember('EX_cpd00528_e',community.mets);
            C = subFlux(loc)/V;
            Vmax_N2l = 0;
            Vmax_N2u = 10000;
            community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00528_e',[Vmax_N2l],'l');
            community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00528_e',[Vmax_N2u],'u');
            % NH2OH
            [~, loc] = ismember('EX_cpd00165_e',community.mets);
            C = subFlux(loc)/V;
            Vmax_NH2OHl = -1.125*C/(C+0.000365); 
            Vmax_NH2OHu = 10000;
            community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00165_e',[Vmax_NH2OHl],'l');
            community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00165_e',[Vmax_NH2OHu],'u');
            
        case 'Nitrobacter winogradskyi'			
			% ammonia
            [TF, loc] = ismember('EX_cpd00013_e',community.mets);
            C = subFlux(loc)/V;
            Vmax_amm = -22.5*C/(C + 0.000464);  
			community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00013_e',[Vmax_amm],'l');		
			% nitrite
            [~, loc] = ismember('EX_cpd00075_e',community.mets);
            C = subFlux(loc)/V;
            Vmax_NO2 = -22.5*C/(C + 0.000464);
            community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00075_e',[Vmax_NO2],'l');
            % NO3
            [~, loc] = ismember('EX_cpd00209_e',community.mets);
            C = subFlux(loc)/V;
            Vmax_NO3 = 0;
			community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00209_e',[Vmax_NO3],'l');            
            % NO
            [~, loc] = ismember('EX_cpd00418_e',community.mets);
            C = subFlux(loc)/V;
            Vmax_NOl = -22.5*C/(C + 0.0518);
            Vmax_NOu=10000;
			community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00418_e',[Vmax_NOl],'l');
			community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00418_e',[Vmax_NOu],'u');      
            % CO2
            [TF, loc] = ismember('EX_cpd00011_e',community.mets);
            C = subFlux(loc)/V;
            Vmax_CO2l = -22.5*C/(C + 0.00037);  
            Vmax_CO2u = 10000;
			community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00011_e',[Vmax_CO2l],'l');
			community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00011_e',[Vmax_CO2u],'u');			
            % bicarb
            [~, loc] = ismember('EX_cpd00242_e',community.mets);
            C = subFlux(loc)/V;
            Vmax_HCO3l = -22.5*C/(C + 0.00037); 
            Vmax_HCO3u=  10000;
            community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00242_e',[Vmax_HCO3l],'l');
            community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00242_e',[Vmax_HCO3u],'u');      
            % O2
            [~, loc] = ismember('EX_cpd00007_e',community.mets);
            C = subFlux(loc)/V
            Vmax_O2l =  -22.5*C/(C + 0.0518);  
            Vmax_O2u = 10000;
			community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00007_e',[Vmax_O2l],'l');
            community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00007_e',[Vmax_O2u],'u');
              % CHO2
            [~, loc] = ismember('EX_cpd00047_e',community.mets);
            C = subFlux(loc)/V;
            Vmax_CHO2l = -11.25*C/(C+0.00037);  
            Vmax_CHO2u = 10000;
            community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00047_e',[Vmax_CHO2l],'l');
            community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00047_e',[Vmax_CHO2u],'u');
            % N2O
            [~, loc] = ismember('EX_cpd00659_e',community.mets);
            C = subFlux(loc)/V;
            Vmax_N2Ol = 0;  
            Vmax_N2Ou = 10000;
			community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00659_e',[Vmax_N2Ol],'l');
            community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00659_e',[Vmax_N2Ou],'u');
            % N2
            [~, loc] = ismember('EX_cpd00528_e',community.mets);
            C = subFlux(loc)/V;
            Vmax_N2l = 0;  
            Vmax_N2u = 0;
            community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00528_e',[Vmax_N2l],'l');
            community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd00528_e',[Vmax_N2u],'u');
        
    end
    
    end