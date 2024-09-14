(* ::Package:: *)

ProjectionOperatorN::usage = "ProjectionOperatorN[group,r,p,s1,ms1,p1,s2,ms2,p2] generates two-hadron operators by projection.";


Begin["`Projection`"];


(* Projection method for hadron in any groups *)
ProjectionOperatorN[ptot_,rep_,r_,momTuple_,parTuple_,msTuple_]:=Module[{group,repO,repparity,ope,Npar,sTuple,pTuple,singleMul,singleOpe,HadronName},
group=MomToGroup[ptot];
(* Separate the parity of the rep. *)
If[MemberQ[{"Oh","OhD"},group],
repO=StringTake[rep,{1,-2}];
repparity=If[StringPart[rep,-1]==="+",1,If[StringPart[rep,-1]==="-",-1,Print["wrong rep. parity"]]],
repO=rep;
repparity=1];
(* Number of particles *)
Npar=Length[momTuple];
(* Extract the particle information *)
sTuple=ConstantArray[0,Npar];
pTuple=ConstantArray[0,Npar];
Do[{sTuple[[i]],pTuple[[i]]}=parTuple[[i]],{i,Npar}];
ope=0;
(* The rotation elements *)
Do[
singleMul=1;
Do[
singleOpe=0;
Do[singleOpe+=H[pari,sTuple[[pari]],msp,pTuple[[pari]]][RotationMatrix[OhD["\[Omega]"][[i]],OhD["n"][[i]]] . momTuple[[pari]]] Conjugate[MyDD[sTuple[[pari]],msp,msTuple[[pari]],i]],{msp,-sTuple[[pari]],sTuple[[pari]]}];
singleMul=singleMul singleOpe,{pari,Npar}];
ope+=Representation[group,repO,0,i,r,r] singleMul,{i,Gele[group]["+"]}];
(* The rest of the little group *)
Do[
singleMul=1;
Do[
singleOpe=0;
Do[singleOpe+=parparity[pTuple[[pari]]] H[pari,sTuple[[pari]],msp,pTuple[[pari]]][-RotationMatrix[OhD["\[Omega]"][[i]],OhD["n"][[i]]] . momTuple[[pari]]] Conjugate[MyDD[sTuple[[pari]],msp,msTuple[[pari]],i]],{msp,-sTuple[[pari]],sTuple[[pari]]}];
singleMul=singleMul singleOpe,{pari,Npar}];
ope+=repparity Representation[group,repO,1,i,r,r] singleMul,{i,Gele[group]["-"]}];
ope=OSimplify[ope];
Return[ope];
];


ProjectionOperatorNString[ptot_,rep_,r_,momTuple_,parTuple_,msTuple_]:=Module[{group,repO,repparity,ope,Npar,sTuple,pTuple,singleMul,singleOpe,HadronName},
group=MomToGroup[ptot];
(* Separate the parity of the rep. *)
If[MemberQ[{"Oh","OhD"},group],
repO=StringTake[rep,{1,-2}];
repparity=If[StringPart[rep,-1]==="+",1,If[StringPart[rep,-1]==="-",-1,Print["wrong rep. parity"]]],
repO=rep;
repparity=1];
(* Number of particles *)
Npar=Length[momTuple];
(* Extract the particle information *)
sTuple=ConstantArray[0,Npar];
pTuple=ConstantArray[0,Npar];
Do[{sTuple[[i]],pTuple[[i]]}=parTuple[[i]],{i,Npar}];
ope=0;
(* The rotation elements *)
Do[
singleMul=1;
Do[
singleOpe=0;
HadronName="H"<>ToString[pari];
Do[singleOpe+=\!\(\*SubscriptBox[\(HadronName\), \({sTuple[\([pari]\)], pTuple[\([pari]\)]}\)]\)[msp][RotationMatrix[OhD["\[Omega]"][[i]],OhD["n"][[i]]] . momTuple[[pari]]] Conjugate[MyDD[sTuple[[pari]],msp,msTuple[[pari]],i]],{msp,-sTuple[[pari]],sTuple[[pari]]}];
singleMul=singleMul singleOpe,{pari,Npar}];
ope+=Representation[group,repO,0,i,r,r] singleMul,{i,Gele[group]["+"]}];
(* The rest of the little group *)
Do[
singleMul=1;
Do[
singleOpe=0;
HadronName="H"<>ToString[pari];
Do[singleOpe+=parparity[pTuple[[pari]]] \!\(\*SubscriptBox[\(HadronName\), \({sTuple[\([pari]\)], pTuple[\([pari]\)]}\)]\)[msp][-RotationMatrix[OhD["\[Omega]"][[i]],OhD["n"][[i]]] . momTuple[[pari]]] Conjugate[MyDD[sTuple[[pari]],msp,msTuple[[pari]],i]],{msp,-sTuple[[pari]],sTuple[[pari]]}];
singleMul=singleMul singleOpe,{pari,Npar}];
ope+=repparity Representation[group,repO,1,i,r,r] singleMul,{i,Gele[group]["-"]}];
ope=OSimplify[ope];
Return[ope];
];


End[];
