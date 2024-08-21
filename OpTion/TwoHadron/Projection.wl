(* ::Package:: *)

ProjectionOperator::usage = "ProjectionOperator[rep,r,p,s1,ms1,p1,s2,ms2,p2] generates two-hadron operators in the O group by projection.";
ProjectionOperatorD::usage = "ProjectionOperatorD[rep,r,p,s1,ms1,p1,s2,ms2,p2] generates two-hadron operators in the \!\(\*SuperscriptBox[\(O\), \(D\)]\) group by projection.";
ProjectionOperatorLittle::usage = "ProjectionOperatorLittle[group,rep,r,mom1,mom2,s1,ms1,p1,s2,ms2,p2] generates two-hadron operators in the little groups by projection.";
ProjectionOperatorLittleD::usage = "ProjectionOperatorLittleD[group,rep,r,mom1,mom2,s1,ms1,p1,s2,ms2,p2] generates two-hadron operators in the double cover of the little groups by projection.";


Begin["`Projection`"];


(* Projection method *)
ProjectionOperator[rep_,r_,p_,s1_,ms1_,p1_,s2_,ms2_,p2_]:=Module[{ope,singleOpe1,singleOpe2,repparity,repO},
repO=StringTake[rep,{1,-2}];
repparity=If[StringPart[rep,-1]==="+",1,If[StringPart[rep,-1]==="-",-1,Print["wrong rep. parity"]]];
ope=0;
(* The rotation elements *)
Do[
singleOpe1=0;
singleOpe2=0;
Do[singleOpe1+=H[1,s1,ms1p,p1][RotationMatrix[Oh["\[Omega]"][[i]],Oh["n"][[i]]] . p] Conjugate[MyD[s1,ms1p,ms1,i]],{ms1p,-s1,s1}];
Do[singleOpe2+=H[2,s2,ms2p,p2][-RotationMatrix[Oh["\[Omega]"][[i]],Oh["n"][[i]]] . p] Conjugate[MyD[s2,ms2p,ms2,i]],{ms2p,-s2,s2}];
ope+=Orep[repO,i,r,r]singleOpe1 singleOpe2 ,{i,1,24}];
(* The rest of the Subscript[O, h] group *)
Do[
singleOpe1=0;
singleOpe2=0;
Do[singleOpe1+=parparity[p1]H[1,s1,ms1p,p1][-RotationMatrix[Oh["\[Omega]"][[i]],Oh["n"][[i]]] . p] Conjugate[MyD[s1,ms1p,ms1,i]],{ms1p,-s1,s1}];
Do[singleOpe2+=parparity[p2]H[2,s2,ms2p,p2][RotationMatrix[Oh["\[Omega]"][[i]],Oh["n"][[i]]] . p] Conjugate[MyD[s2,ms2p,ms2,i]],{ms2p,-s2,s2}];
ope+=repparity Orep[repO,i,r,r]singleOpe1 singleOpe2 ,{i,1,24}];
ope=OSimplify[ope];
Return[ope];
];


(* Projection method *)
ProjectionOperatorD[rep_,r_,p_,s1_,ms1_,p1_,s2_,ms2_,p2_]:=Module[{ope,singleOpe1,singleOpe2,repparity,repO},
repO=StringTake[rep,{1,-2}];
repparity=If[StringPart[rep,-1]==="+",1,If[StringPart[rep,-1]==="-",-1,Print["wrong rep. parity"]]];
ope=0;
(* The rotation elements *)
Do[
singleOpe1=0;
singleOpe2=0;
Do[singleOpe1+=H[1,s1,ms1p,p1][RotationMatrix[OhD["\[Omega]"][[i]],OhD["n"][[i]]] . p] Conjugate[MyDD[s1,ms1p,ms1,i]],{ms1p,-s1,s1}];
Do[singleOpe2+=H[2,s2,ms2p,p2][-RotationMatrix[OhD["\[Omega]"][[i]],OhD["n"][[i]]] . p] Conjugate[MyDD[s2,ms2p,ms2,i]],{ms2p,-s2,s2}];
ope+=ODrep[repO,i,r,r]singleOpe1 singleOpe2 ,{i,1,48}];
(* The rest of the Subscript[O, h] group *)
Do[
singleOpe1=0;
singleOpe2=0;
Do[singleOpe1+=parparity[p1]H[1,s1,ms1p,p1][-RotationMatrix[OhD["\[Omega]"][[i]],OhD["n"][[i]]] . p] Conjugate[MyDD[s1,ms1p,ms1,i]],{ms1p,-s1,s1}];
Do[singleOpe2+=parparity[p2]H[2,s2,ms2p,p2][RotationMatrix[OhD["\[Omega]"][[i]],OhD["n"][[i]]] . p] Conjugate[MyDD[s2,ms2p,ms2,i]],{ms2p,-s2,s2}];
ope+=repparity ODrep[repO,i,r,r]singleOpe1 singleOpe2 ,{i,1,48}];
ope=OSimplify[ope];
Return[ope];
];


(* Projection method for little groups *)
ProjectionOperatorLittle[group_,rep_,r_,mom1_,mom2_,s1_,ms1_,p1_,s2_,ms2_,p2_]:=Module[{ope,singleOpe1,singleOpe2},
(* No rep. parity *)
ope=0;
(* The rotation elements *)
Do[
singleOpe1=0;
singleOpe2=0;
Do[singleOpe1+=H[1,s1,ms1p,p1][RotationMatrix[Oh["\[Omega]"][[i]],Oh["n"][[i]]] . mom1] Conjugate[MyD[s1,ms1p,ms1,i]],{ms1p,-s1,s1}];
Do[singleOpe2+=H[2,s2,ms2p,p2][RotationMatrix[Oh["\[Omega]"][[i]],Oh["n"][[i]]] . mom2] Conjugate[MyD[s2,ms2p,ms2,i]],{ms2p,-s2,s2}];
ope+=Littlerep[group,rep,0,i,r,r]singleOpe1 singleOpe2 ,{i,Gele[group]["+"]}];
(* The rest of the little group *)
Do[
singleOpe1=0;
singleOpe2=0;
Do[singleOpe1+=parparity[p1]H[1,s1,ms1p,p1][-RotationMatrix[Oh["\[Omega]"][[i]],Oh["n"][[i]]] . mom1] Conjugate[MyD[s1,ms1p,ms1,i]],{ms1p,-s1,s1}];
Do[singleOpe2+=parparity[p2]H[2,s2,ms2p,p2][-RotationMatrix[Oh["\[Omega]"][[i]],Oh["n"][[i]]] . mom2] Conjugate[MyD[s2,ms2p,ms2,i]],{ms2p,-s2,s2}];
(* Subscript[\[Sigma], x]=Subscript[IC, x](\[Pi]) *)
(* After Subscript[C, x](\[Pi]) or Subscript[C, x,-y](\[Pi])*)
ope+=Littlerep[group,rep,1,i,r,r]singleOpe1 singleOpe2 ,{i,Gele[group]["-"]}];
ope=OSimplify[ope];
Return[ope];
];


(* Projection method for little groups *)
ProjectionOperatorLittleD[group_,rep_,r_,mom1_,mom2_,s1_,ms1_,p1_,s2_,ms2_,p2_]:=Module[{ope,singleOpe1,singleOpe2},
(* No rep. parity *)
ope=0;
(* The rotation elements *)
Do[
singleOpe1=0;
singleOpe2=0;
Do[singleOpe1+=H[1,s1,ms1p,p1][RotationMatrix[OhD["\[Omega]"][[i]],OhD["n"][[i]]] . mom1] Conjugate[MyDD[s1,ms1p,ms1,i]],{ms1p,-s1,s1}];
Do[singleOpe2+=H[2,s2,ms2p,p2][RotationMatrix[OhD["\[Omega]"][[i]],OhD["n"][[i]]] . mom2] Conjugate[MyDD[s2,ms2p,ms2,i]],{ms2p,-s2,s2}];
ope+=LittlerepD[group,rep,0,i,r,r]singleOpe1 singleOpe2 ,{i,Gele[group]["+"]}];
(* The rest of the little group *)
Do[
singleOpe1=0;
singleOpe2=0;
Do[singleOpe1+=parparity[p1]H[1,s1,ms1p,p1][-RotationMatrix[OhD["\[Omega]"][[i]],OhD["n"][[i]]] . mom1] Conjugate[MyDD[s1,ms1p,ms1,i]],{ms1p,-s1,s1}];
Do[singleOpe2+=parparity[p2]H[2,s2,ms2p,p2][-RotationMatrix[OhD["\[Omega]"][[i]],OhD["n"][[i]]] . mom2] Conjugate[MyDD[s2,ms2p,ms2,i]],{ms2p,-s2,s2}];
(* Subscript[\[Sigma], x]=Subscript[IC, x](\[Pi]) *)
(* After Subscript[C, x](\[Pi]) or Subscript[C, x,-y](\[Pi])*)
ope+=LittlerepD[group,rep,1,i,r,r]singleOpe1 singleOpe2 ,{i,Gele[group]["-"]}];
ope=OSimplify[ope];
Return[ope];
];


End[];
