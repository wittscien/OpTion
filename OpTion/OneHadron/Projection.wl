(* ::Package:: *)

ProjectionOperator1::usage = "ProjectionOperator1[ptot,rep,r,J,M,JD,J13,TG,ND] generates one-hadron operators by projection.";


Begin["`Projection`"];


(* Projection method for hadron in any groups *)
ProjectionOperator1[ptot_,rep_,r_,J_,M_,JD_,J13_,TG_,ND_]:=Module[{group,repO,repparity,P,PGD,ope,singleOpe},
group=MomToGroup[ptot];
(* Separate the parity of the rep. *)
If[MemberQ[{"Oh","OhD"},group],
repO=StringTake[rep,{1,-2}];
repparity=If[StringPart[rep,-1]==="+",1,If[StringPart[rep,-1]==="-",-1,Print["wrong rep. parity"]]],
repO=rep;
repparity=1];
P=If[TG==="S" || TG==="A",1,If[TG==="P" || TG==="V",-1,0]];
PGD=P (-1)^ND;
(* Extract the particle information *)
ope=0;
(* The rotation elements *)
Do[
singleOpe=0;
Do[singleOpe+=OJM[J,msp,JD,J13,TG,ND]Conjugate[MyDD[J,msp,M,i]],{msp,-J,J}];
ope+=Representation[group,repO,0,i,r,r]singleOpe,{i,Gele[group]["+"]}];
(* The rest of the little group *)
Do[
singleOpe=0;
Do[singleOpe+=PGD OJM[J,msp,JD,J13,TG,ND]Conjugate[MyDD[J,msp,M,i]],{msp,-J,J}];
ope+=repparity Representation[group,repO,1,i,r,r]singleOpe,{i,Gele[group]["-"]}];
ope=NCDSimplify[OSimplify[ope]];
Return[ope];
];


End[];
