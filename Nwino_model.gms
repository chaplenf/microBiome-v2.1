$ontext
LP model to be called by tr3_run.m.

In this model we assume the data for the LP is passed in via a call to the
gams() Mex-function.
$offtext

$set matout "'A:\Neuro_Supplementary_Materials\matsol.gdx', returnStat, Z, v "

Sets
Reactions
/
$include A:\Neuro_Supplementary_Materials\Nwi_Reactions.txt
/
Metabolites
/
$include A:\Neuro_Supplementary_Materials\Nwi_Compounds.txt
/
;
Parameter
Vmax(Reactions)
/
$include A:\Neuro_Supplementary_Materials\Nwi_UpperBounds.txt
/
;


Parameter
Vmin(Reactions)
/
$Include A:\Neuro_Supplementary_Materials\Nwi_LowerBounds.txt
/
;

$include A:\Neuro_Supplementary_Materials\Nwi_Objective.gms

$include A:\Neuro_Supplementary_Materials\iFC579.gms

Variables
v(Reactions)
Z;

Equations
massbalance(Metabolites) S.v = 0
objectivefunction        Z=c.v

NitrogenBalance;

massbalance(Metabolites).. sum(Reactions, S(Metabolites,Reactions)*v(Reactions))=e=0;
objectivefunction .. Z=e=sum(Reactions,objfunc(Reactions)*v(Reactions));
NitrogenBalance .. 2*v('R1094')+2*v('R1095')+v('R1078')+v('R1088')=l=-(v('R1076'));


*Place lower bounds of 0 for all fluxes
v.lo(Reactions)=Vmin(Reactions);
*Place upper bounds of 10 for all fluxes
v.up(Reactions)=Vmax(Reactions);

model FBA /all/;
solve FBA using lp maximizing Z;

set stat /modelstat,solvestat/;
parameter returnStat(stat);
returnStat('modelstat') = FBA.modelstat;
returnStat('solvestat') = FBA.solvestat;

execute_unload %matout%;
