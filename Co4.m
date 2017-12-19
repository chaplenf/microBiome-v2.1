function [t, V, Vg, X, L, LNM, G, AB, BB, R] = Co4(community, bucket, tspan, initialConditions)

%   tspan: simulation time
%   intialConditions: the initial reactor volume, biomass and metabolite concentrations

    numSpecies = length(community.species);
    numLiquidSpecies = length(community.Lfeed);
    numLNM = length(community.LNMfeed);
    numGasSpecies = length(community.Gfeed);
    numABucket = length(bucket.numabucket);
    numBBucket=length(bucket.numbbucket);
    numR=length(bucket.flux);


        DO=@coCulture_ode4;
        options=odeset('NonNegative',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 17 18 19 20 21 22 23 24 25 26]);
        [t,y] = ode15s(DO,tspan,initialConditions,options);


V = y(:,1);
Vg = y(:,2);
X = y(:,3:2+numSpecies);
L = y(:,2+numSpecies+1: 2+numSpecies+numLiquidSpecies);
LNM = y(:,2+numSpecies+numLiquidSpecies+1: 2+numSpecies+numLiquidSpecies+numLNM);
G = y(:,2+numSpecies+numLiquidSpecies+numLNM+1:2+numSpecies+numLiquidSpecies+numLNM+numGasSpecies);
AB = y(:,2+numSpecies+numLiquidSpecies+numLNM+numGasSpecies+1:2+numSpecies+numLiquidSpecies+numLNM+numGasSpecies+numABucket);
BB = y(:,2+numSpecies+numLiquidSpecies+numLNM+numGasSpecies+numABucket+1:2+numSpecies+numLiquidSpecies+numLNM+numGasSpecies+numABucket+numBBucket);
R = y(:,2+numSpecies+numLiquidSpecies+numLNM+numGasSpecies+numABucket+numBBucket+1:2+numSpecies+numLiquidSpecies+numLNM+numGasSpecies+numABucket+numBBucket+numR);        
end

