(* ::Package:: *)

(*OJM::usage = "OJM[J,M,JD,J13,TG,ND] constructs |J,M> one-hadron operators in the rest frame and infinite volume.";*)
(*OJ\[Lambda]::usage = "OJ\[Lambda][group,J, \[Lambda],JD,J13,TG,ND] constructs |J,\[Lambda]> one-hadron operators in the moving frame from |J,M> in infinite volum.";*)
OSub::usage = "OSub[group_,rep_,r_,J_, \[Lambda]_,JD_,J13_,TG_,ND_] constructs one-hadron operators for the little groups.";


Begin["`Projection`"];


(* ND settings *)
(* Note on 2023.11.28: I deleted most of the code in the NC package and commented all warnings to make it lite. *)
CommuteEverything[];
SetNonCommutative[D];


(* O^JM *)
(* Note that here we use the subscript notation but function notation for particles. *)
OJM[J_,M_,JD_,J13_,TG_,ND_]:=Module[{ope,m2,m3},
(* TG: type of \[CapitalGamma]; ND: number of D *)
ope=0;
If[TG==="V" || TG==="A",
If[ND===0,If[J===1,ope=Subscript[V, M]]];
If[ND===1,Do[m2=M-m1;If[-1<=m2<=1,ope+=MyClebschGordan[{1,m1},{1,m2},{J,M}]Subscript[V, m1] Subscript[D, m2]],{m1,-1,1}]];
If[ND===2,Do[Do[m3=M-mD;m2=mD-m1;If[-1<=m2<=1 && -1<=m3<=1,ope+=MyClebschGordan[{1,m3},{JD,mD},{J,M}]MyClebschGordan[{1,m1},{1,m2},{JD,mD}]Subscript[V, m3]Subscript[D, m1]**Subscript[D, m2]],{m1,-1,1}],{mD,-JD,JD}]];
];
If[TG==="S" || TG==="P",
If[ND===0,If[J===0,ope=S]];
If[ND===1,If[J===1,ope=S Subscript[D, M]]];
If[ND===2,Do[m2=M-m1;If[-1<=m2<=1,ope+=MyClebschGordan[{1,m1},{1,m2},{J,M}]S Subscript[D, m1]**Subscript[D, m2]],{m1,-1,1}]];
];
If[TG==="A",ope=ope//.{V->A}];
If[TG==="P",ope=ope//.{S->P}];
Return[ope];
];


(* O^J\[Lambda] *)
OJ\[Lambda][group_,J_, \[Lambda]_,JD_,J13_,TG_,ND_]:=Module[{ope,\[Psi],\[Theta],\[Phi]},
ope=0;
(* \[Psi],\[Theta],\[Phi] have taken a minus sign *)
If[group==="C4v",\[Psi]=0;\[Theta]=0;\[Phi]=0,If[group==="C2v",\[Psi]=-(\[Pi]/2);\[Theta]=-(\[Pi]/4);\[Phi]=\[Pi]/2,If[group==="C3v",\[Psi]=-(\[Pi]/4);\[Theta]=-ArcCos[1/Sqrt[3]];\[Phi]=0]]];
Do[ope+=Conjugate[WignerD[{J,M,\[Lambda]},\[Psi],\[Theta],\[Phi]]]OJM[J,M,JD,J13,TG,ND],{M,-J,J}];
Return[ope];
];


(* Rest frame: Subsuperscript[O, \[CapitalLambda],\[Mu], [J]] *)
OSub[rep_,r_,J_,JD_,J13_,TG_,ND_]:=Module[{ope,repparity,repO,P,PGD},
repO=StringTake[rep,{1,-2}];
repparity=If[StringPart[rep,-1]==="+",1,If[StringPart[rep,-1]==="-",-1,Print["wrong rep. parity"]]];
P=If[TG==="S" || TG==="A",1,If[TG==="P" || TG==="V",-1,0]];
PGD=P (-1)^ND; (* Parity of Gamma\[Cross]D *)
If[repparity=!=PGD,Print["Wrong parity"];Return[0]];
ope=0;
Do[ope+=subduction[{J,repO}][[r,-M+J+1]]OJM[J,M,JD,J13,TG,ND],{M,-J,J}];
ope=NCDSimplify[ope];
Return[ope];
];


(* Moving frame: Subsuperscript[O, \[CapitalLambda],\[Mu], [J,|\[Lambda]|]] *)
OSub[group_,rep_,r_,J_, \[Lambda]_,JD_,J13_,TG_,ND_]:=Module[{ope,P,\[Eta]tilde,PGD},
If[J<Abs[\[Lambda]],Print["Wrong: J < \[Lambda]"];Return[0]];
P=If[TG==="S" || TG==="A",1,If[TG==="P" || TG==="V",-1,0]];
PGD=P (-1)^ND; (* Parity of Gamma\[Cross]D *)
\[Eta]tilde=PGD (-1)^J;
ope=0;
Do[ope+=subductionLittle[group,rep,r,\[Eta]tilde,\[Lambda]hat]OJ\[Lambda][group,J,\[Lambda]hat,JD,J13,TG,ND],{\[Lambda]hat,{Abs[\[Lambda]],-Abs[\[Lambda]]}}];
ope=NCDSimplify[ope];
Return[ope];
];


End[];
