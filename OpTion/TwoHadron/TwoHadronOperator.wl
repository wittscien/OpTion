(* ::Package:: *)

TwoHadronOperator::usage = "TwoHadronOperator[ptot,rep,r,mom,par1,ms1,par2,ms2] generates two-hadron operators by projection.";
TwoHadronOperatorPartialWave::usage = "TwoHadronOperatorPartialWave[ptot,rep,r,mom,J,L,S,par1,par2] generates two-hadron operators by partial wave coupling."
TwoHadronOperatorAll::usage = "TwoHadronOperatorAll[ptot,rep,r,maxmom,par1,par2] generates all possible two-hadron operators by projection.";
TwoHadronOperatorPartialWaveAll::usage = "TwoHadronOperatorPartialWaveAll[ptot,rep,r,maxmom,par1,par2,L] generates all possible two-hadron operators by partial wave coupling.";


Begin["`TwoHadronOperator`"];


(* Combine all functions *)
TwoHadronOperator[ptot_,rep_,r_,mom_,par1_,par2_,ms1_,ms2_]:=Module[{s1,p1,s2,p2,group,op},
(* Translates the total momentum into the (little) group *)
(* Add arbitrary spin *)
If[Head[par1]===List,{s1,p1}=par1,{s1,p1}=ParToSpinParity[par1]];
If[Head[par2]===List,{s2,p2}=par2,{s2,p2}=ParToSpinParity[par2]];
group=If[IntegerQ[s1+s2],MomToGroup[ptot],MomToGroup[ptot,1]];
(* Distribute the funtions *)
If[IntegerQ[s1+s2] && ptot === {0,0,0},op=ProjectionOperator[rep,r,mom,s1,ms1,p1,s2,ms2,p2]];
If[!IntegerQ[s1+s2] && ptot === {0,0,0},op=ProjectionOperatorD[rep,r,mom,s1,ms1,p1,s2,ms2,p2]];
If[IntegerQ[s1+s2] && ptot =!= {0,0,0},op=ProjectionOperatorLittle[group,rep,r,mom,ptot-mom,s1,ms1,p1,s2,ms2,p2]];
If[!IntegerQ[s1+s2] && ptot =!= {0,0,0},op=ProjectionOperatorLittleD[group,rep,r,mom,ptot-mom,s1,ms1,p1,s2,ms2,p2]];
Return[op];
];


(* Partial-wave method *)
TwoHadronOperatorPartialWave[ptot_,rep_,r_,mom_,J_,L_,S_,par1_,par2_]:=Module[{s1,p1,s2,p2,group,op},
(* Translates the total momentum into the (little) group *)
(* Add arbitrary spin *)
If[Head[par1]===List,{s1,p1}=par1,{s1,p1}=ParToSpinParity[par1]];
If[Head[par2]===List,{s2,p2}=par2,{s2,p2}=ParToSpinParity[par2]];
group="Oh";
(* Distribute the funtions *)
If[ptot === {0,0,0},op=PartialWaveOperatorProjected[rep,r,mom,J,L,S,s1,p1,s2,p2],Throw["Only for rest frame"]];
Return[op];
];


(* Search for all operators *)
TwoHadronOperatorAll[ptot_,rep_,r_,maxmom_,par1_,par2_]:=Module[{s1,p1,s2,p2,group,momList,opList,op},
(* Translates the total momentum into the (little) group *)
If[Head[par1]===List,{s1,p1}=par1,{s1,p1}=ParToSpinParity[par1]];
If[Head[par2]===List,{s2,p2}=par2,{s2,p2}=ParToSpinParity[par2]];
group=If[IntegerQ[s1+s2],MomToGroup[ptot],MomToGroup[ptot,1]];
momList=SortBy[Select[GenerateMomVectors[maxmom],Norm[ptot - #] <= maxmom &],Norm[#]^2 + Norm[ptot - #]^2 &];
(* Search for all possible spins and momenta *)
opList={};
Do[
op=TwoHadronOperator[ptot,rep,r,mom,par1,par2,ms1,ms2];
If[!MemberQ[opList,op] && op =!= 0 && IsLinearlyIndependent[opList,op],AppendTo[opList,op]],{mom,momList},{ms1,-s1,s1},{ms2,-s2,s2}];
Return[opList];
];

TwoHadronOperatorAllOld[ptot_,rep_,r_,maxmom_,par1_,par2_]:=Module[{s1,p1,s2,p2,group,momList,opList,op},
(* Translates the total momentum into the (little) group *)
(* Add arbitrary spin *)
If[Head[par1]===List,{s1,p1}=par1,{s1,p1}=ParToSpinParity[par1]];
If[Head[par2]===List,{s2,p2}=par2,{s2,p2}=ParToSpinParity[par2]];
group=If[IntegerQ[s1+s2],MomToGroup[ptot],MomToGroup[ptot,1]];
momList=SortBy[Select[GenerateMomVectors[maxmom],Norm[ptot - #] <= maxmom &],Norm[#]^2 + Norm[ptot - #]^2 &];
(* Search for all possible spins and momenta *)
opList={};
Do[
If[IntegerQ[s1+s2] && ptot === {0,0,0},op=ProjectionOperator[rep,r,mom,s1,ms1,p1,s2,ms2,p2]];
If[!IntegerQ[s1+s2] && ptot === {0,0,0},op=ProjectionOperatorD[rep,r,mom,s1,ms1,p1,s2,ms2,p2]];
If[IntegerQ[s1+s2] && ptot =!= {0,0,0},op=ProjectionOperatorLittle[group,rep,r,mom,ptot-mom,s1,ms1,p1,s2,ms2,p2]];
If[!IntegerQ[s1+s2] && ptot =!= {0,0,0},op=ProjectionOperatorLittleD[group,rep,r,mom,ptot-mom,s1,ms1,p1,s2,ms2,p2]];
If[!MemberQ[opList,op] && op =!= 0,AppendTo[opList,op]],{mom,momList},{ms1,-s1,s1},{ms2,-s2,s2}];
Return[opList];
];


(* Search for all operators *)
TwoHadronOperatorPartialWaveAll[ptot_,rep_,r_,maxmom_,par1_,par2_,L_]:=Module[{repO,repparity,s1,p1,s2,p2,group,momList,opList,op},
(* Translates the total momentum into the (little) group *)
If[Head[par1]===List,{s1,p1}=par1,{s1,p1}=ParToSpinParity[par1]];
If[Head[par2]===List,{s2,p2}=par2,{s2,p2}=ParToSpinParity[par2]];
repO=StringTake[rep,{1,-2}];
repparity=If[StringPart[rep,-1]==="+",1,If[StringPart[rep,-1]==="-",-1,Print["wrong rep. parity"]]];
If[repparity=!=(-1)^L parparity[p1] parparity[p2],Print["The choice of L does not match the parity of the representation. What makes you think so?"];Return[{}]];
group="Oh";
momList=SortBy[Select[GenerateMomVectors[maxmom],Norm[ptot - #] <= maxmom &],Norm[#]^2 + Norm[ptot - #]^2 &];
(* Search for all possible spins and momenta *)
opList={};
Do[
op=TwoHadronOperatorPartialWave[ptot,repO,r,mom,J,L,S,par1,par2];
If[!MemberQ[opList,op] && op =!= 0,AppendTo[opList,op]],{mom,momList},{S,Abs[s1-s2],s1+s2},{J,Abs[L-S],L+S}];
Return[opList];
];


TwoHadronOperatorAll[rep_,r_,maxmom_,par1_,par2_]:=TwoHadronOperatorAll[{0,0,0},rep,r,maxmom,par1,par2];


End[];
