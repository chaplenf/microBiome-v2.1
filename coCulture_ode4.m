function dy = coCulture_ode4(t, y, ~)
global community bucket
%ODE framework
    
dy = zeros(size(y));
    numSpecies = length(community.species);
    numLiquidSpecies = length(community.Lfeed);
    numLiquidNMSpecies = length(community.LNMfeed);
    numGasSpecies = length(community.Gfeed);
    numLB = length(bucket.liquidreactormassbalance);
    numLNMB=length(bucket.NMliquidreactormassbalance);
    numGB = length(bucket.gasreactormassbalance);
    numMT = length(bucket.masstransfer);
    numBBuckets = length(bucket.numbbucket);
    numABuckets = length(bucket.numabucket);
    numInternalFluxes = length(bucket.flux);
V = y(1);                     %Bioreactor liquid volume [L]
Vg = y(2);                    %Bioreactor gas volume [L]
for j = 1:numLiquidSpecies
    L(j) = y(2+numSpecies+j);  %Substrates [mmol]
end
% assigning growth rates and metabolic production/consumption rates
    % here, the rates are calculated using FBA through linear programming models 
    % run in the GAMS environment

vs =zeros(numSpecies,((numInternalFluxes+numBBuckets)/numSpecies));
mu = zeros(1, numSpecies);

reactionIdentifiers = [community.mets bucket.reac];
subFlux = [L];

uptakeRates4;  % updating FBA organism uptake constraints

% calculating the growth rates and production/consumption rates of the organisms

for i = 1:1
%try    
switch community.species{1,i}.description
        case 'Nitrosomonas europaea'
        Reac=textread('Neu_Reactions.txt','%s');
        Reacsize=size(Reac);
        filename = 'Neu_UpperBounds.txt';
        fid = fopen(filename, 'w');
        for row=1:Reacsize
        fprintf(fid,'R%d %f\n',row,community.species{1,i}.ub(row,:));        
        end
        fclose(fid);
        filename = 'Neu_LowerBounds.txt';
        fid = fopen(filename, 'w');
        for row=1:Reacsize
        fprintf(fid,'R%d %f\n',row,community.species{1,i}.lb(row,:)); 
        end
        fclose(fid);
    [sNeu,~,vNeu] = gams('Neuro_model1.gms');
        r=sNeu.val(1,2);
            if (r ~= 1)    
                mu(1) = 0;
               community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd11416(e)',0,'l');  
               community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd11416(e)',0,'u');             
                disp('deathphase');
              else
                xfake=vNeu.val(:,1);
                yfake=vNeu.val(:,2);
                [~,reacIDModel]=ismember(reactionIdentifiers,community.species{1,i}.rxns);
                sizefake=size(xfake);
                [~,loc]=ismember('EX_cpd11416(e)',community.species{1,i}.rxns);   
                for fake=1:sizefake
                   for j=1:((numInternalFluxes+numBBuckets)/numSpecies)
                        if (xfake(fake)==reacIDModel(j))
                           vs(i,j) = yfake(fake);
                        end  
                        if (loc==xfake(fake))
                           mu(i) = yfake(fake);
                        end     
                    end
                end
            end
       case 'Nitrobacter winogradskyi'
        Reac=textread('Nwi_Reactions.txt','%s');
        Reacsize=size(Reac);
        filename = 'Nwi_UpperBounds.txt';
        fid = fopen(filename, 'w');
        for row=1:Reacsize
        fprintf(fid,'R%d %f\n',row,community.species{1,i}.ub(row,:));        
        end
        fclose(fid);
        filename = 'Nwi_LowerBounds.txt';
        fid = fopen(filename, 'w');
        for row=1:Reacsize
        fprintf(fid,'R%d %f\n',row,community.species{1,i}.lb(row,:)); 
        end
        fclose(fid);
       [sNwi,~,vNwi] = gams('Nwino_model.gms');
        r=sNwi.val(1,2);
            if (r ~= 1)    
                mu(2) = 0;
              community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd11416(e)',0,'l');  
              community.species{1,i} = changeRxnBounds(community.species{1,i},'EX_cpd11416(e)',0,'u');             
                disp('deathphase');
              else
                xNwifake=vNwi.val(:,1);
                yNwifake=vNwi.val(:,2);
                [~,reacIDModel]=ismember(reactionIdentifiers,community.species{1,i}.rxns);
                sizefake=size(xNwifake);
                [~,loc]=ismember('EX_cpd11416(e)',community.species{1,i}.rxns);   
                for fake=1:sizefake
                   for j=1:((numInternalFluxes+numBBuckets)/numSpecies)
                        if (xNwifake(fake)==reacIDModel(j))
                           vs(i,j) = yNwifake(fake);
                        end  
                        if (loc==xNwifake(fake))
                           mu(i) = yNwifake(fake);
                        end     
                    end
                end
            end
end
%catch
%    mu(i) = 0;
%   disp('catch');
%end
end
HONOModelODE;  %updating abiotic and biotic reaction sets
t, vs(:,:)








