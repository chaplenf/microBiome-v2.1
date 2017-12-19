$ontext
LP model to be called by tr3_run.m.

In this model we assume the data for the LP is passed in via a call to the
gams() Mex-function.
$offtext

$set matout "'A:\Neuro_Supplementary_Materials\matsol.gdx', returnStat, Z, v "

Sets
Reactions
/
$include A:\Neuro_Supplementary_Materials\Neu_Reactions.txt
/
;
Sets
Metabolites
/
$include A:\Neuro_Supplementary_Materials\Neu_Compounds.txt
/
;
Parameter
Vmax(Reactions)
/
$include A:\Neuro_Supplementary_Materials\Neu_UpperBounds.txt
/
;


Parameter
Vmin(Reactions)
/
$Include A:\Neuro_Supplementary_Materials\Neu_LowerBounds.txt
/
;

$include A:\Neuro_Supplementary_Materials\Neu_Objective.gms

$include A:\Neuro_Supplementary_Materials\iFC578.gms

Variables
v(Reactions)
Z;

Equations
massbalance(Metabolites) S.v = 0
objectivefunction        Z=c.v
NH2OH
HAONOReleaseToMedia
NitrogenBalance
NIRK
NOR
HAON2O
HAON2
P460A
P460B
;

massbalance(Metabolites).. sum(Reactions, S(Metabolites,Reactions)*v(Reactions))=e=0;
objectivefunction .. Z=e=sum(Reactions,objfunc(Reactions)*v(Reactions));
NH2OH ..  v('R1084')=l=0;
NIRK ..  v('R982')=l=1.3;
NOR .. v('R978')=l=0.4;
HAON2O .. v('R987')=e=0;
HAON2 .. v('R988')=g=0;
P460A .. v('R992')=l=0.05;
P460B .. v('R993')=l=0.05;
NitrogenBalance .. v('R1081')+v('R1082')+2*v('R1083')+v('R1084')+2*v('R1087')+v('R1089')+v('R1063')+v('R1015')=l=0;
HAONOReleaseToMedia .. v('R991')+v('R1082')=l=0.025*v('R985');

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
