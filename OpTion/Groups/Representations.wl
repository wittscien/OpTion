(* ::Package:: *)

Orep::usage = "Orep[rep,i,a,b] gives the representation matrix elements [a,b] of \!\(\*SuperscriptBox[\(i\), \(th\)]\) element of the O group.";
ODrep::usage = "ODrep[rep,i,a,b] gives the representation matrix elements [a,b] of \!\(\*SuperscriptBox[\(i\), \(th\)]\) element of the \!\(\*SuperscriptBox[\(O\), \(D\)]\) group.";
Littlerep::usage = "Littlerep[group,rep,inv,i,a,b] gives the representation matrix elements [a,b] of \!\(\*SuperscriptBox[\(i\), \(th\)]\) element of little groups.";
LittlerepD::usage = "LittlerepD[group,rep,inv,i,a,b] gives the representation matrix elements [a,b] of \!\(\*SuperscriptBox[\(i\), \(th\)]\) element of the double cover of the little groups.";
Representation::usage = "Representation[group,rep,inv,i,a,b] gives the representation matrix elements [a,b] of \!\(\*SuperscriptBox[\(i\), \(th\)]\) element of the group.";


Begin["`Representations`"];


(* O group representation. Reps. like Subsuperscript[A, 1, -] will be handled in ProjectionOperator. *)
Orep[rep_,i_,a_,b_]:=Module[{Id,mat,matele},
matele=0;
Id={{1,0},{0,1}};
If[rep==="A1",matele=1];
If[rep==="T1" || rep==="T2",matele=Cos[Oh["\[Omega]"][[i]]] KroneckerDelta[a,b]+(1-Cos[Oh["\[Omega]"][[i]]] )Oh["n"][[i]][[a]] Oh["n"][[i]][[b]]-Sin[Oh["\[Omega]"][[i]]] (LeviCivitaTensor[3,List] . Oh["n"][[i]])[[a,b]]];
If[rep ==="T2" && 10<=i<=21,matele=-matele];
If[rep==="E",
If[MemberQ[{1,22,23,24},i],mat=Id];
If[MemberQ[{14,15,18,19},i],mat=PauliMatrix[3]];
If[MemberQ[{2,5,6,9},i],mat=-Cos[\[Pi]/3]Id+I Sin[\[Pi]/3] PauliMatrix[2]];
If[MemberQ[{3,4,7,8},i],mat=-Cos[\[Pi]/3]Id-I Sin[\[Pi]/3] PauliMatrix[2]];
If[MemberQ[{10,11,16,17},i],mat=-Cos[\[Pi]/3]PauliMatrix[3]-Sin[\[Pi]/3] PauliMatrix[1]];
If[MemberQ[{12,13,20,21},i],mat=-Cos[\[Pi]/3]PauliMatrix[3]+Sin[\[Pi]/3] PauliMatrix[1]];
matele=mat[[a,b]]];
Return[matele];
];


ODrep[rep_,i_,a_,b_]:=Module[{Id,mat,matele},
matele=0;
Id={{1,0},{0,1}};
If[rep==="A1",matele=1];
If[rep==="T1" || rep==="T2",matele=Cos[OhD["\[Omega]"][[i]]] KroneckerDelta[a,b]+(1-Cos[OhD["\[Omega]"][[i]]] )OhD["n"][[i]][[a]] OhD["n"][[i]][[b]]-Sin[OhD["\[Omega]"][[i]]] (LeviCivitaTensor[3,List] . OhD["n"][[i]])[[a,b]]];
If[rep ==="T2" && (10<=i<=21 || 10+24<=i<=21+24),matele=-matele];
If[rep==="E",
If[MemberQ[Join[{1,22,23,24},{1,22,23,24}+24],i],mat=Id];
If[MemberQ[Join[{14,15,18,19},{14,15,18,19}+24],i],mat=PauliMatrix[3]];
If[MemberQ[Join[{2,5,6,9},{2,5,6,9}+24],i],mat=-Cos[\[Pi]/3]Id+I Sin[\[Pi]/3] PauliMatrix[2]];
If[MemberQ[Join[{3,4,7,8},{3,4,7,8}+24],i],mat=-Cos[\[Pi]/3]Id-I Sin[\[Pi]/3] PauliMatrix[2]];
If[MemberQ[Join[{10,11,16,17},{10,11,16,17}+24],i],mat=-Cos[\[Pi]/3]PauliMatrix[3]-Sin[\[Pi]/3] PauliMatrix[1]];
If[MemberQ[Join[{12,13,20,21},{12,13,20,21}+24],i],mat=-Cos[\[Pi]/3]PauliMatrix[3]+Sin[\[Pi]/3] PauliMatrix[1]];
matele=mat[[a,b]]];
If[rep==="G1" || rep==="G2",matele=Simplify[MatrixExp[-(I/2)(OhD["n"][[i]][[1]] PauliMatrix[1]+OhD["n"][[i]][[2]] PauliMatrix[2]+OhD["n"][[i]][[3]] PauliMatrix[3])OhD["\[Omega]"][[i]]]][[a,b]]];
If[rep ==="G2" && (10<=i<=21 || 34<=i<=45),matele=-matele];
If[rep==="H",matele=Simplify[MatrixExp[-(I/2)(OhD["n"][[i]][[1]] PauliMatrix3d2[1]+OhD["n"][[i]][[2]] PauliMatrix3d2[2]+OhD["n"][[i]][[3]] PauliMatrix3d2[3])OhD["\[Omega]"][[i]]]][[a,b]]];
Return[matele];
];


(* Little group representation. *)
(* [HSC] PRD 2012 Dudek pipi *)
Littlerep[group_,rep_,inv_,i_,a_,b_]:=Module[{mat,matele},
matele=10000;
(* [0,0,n] *)
(* 1:0, 14:-(\[Pi]/2) 15:\[Pi]/2, 24:\[Pi]; C: 22:0, 19:-(\[Pi]/2) 18:\[Pi]/2, 23:\[Pi] *)
If[group==="C4v",
If[rep==="A1",matele=1];
If[rep==="A2",matele=If[inv===0,1,-1]];
If[rep==="E2",If[inv===0,
If[i===1,mat={{1,0},{0,1}}];
If[i===24,mat={{-1,0},{0,-1}}];
If[i===14,mat={{0,I},{I,0}}];
If[i===15,mat={{0,-I},{-I,0}}];
matele=mat[[a,b]],
If[i===22,mat={{1,0},{0,-1}}];
If[i===23,mat={{-1,0},{0,1}}];
If[i===19,mat={{0,-I},{I,0}}];
If[i===18,mat={{0,I},{-I,0}}];
matele=mat[[a,b]]]];
If[rep==="B1",matele=If[i===1 || i===24 || i===22 || i===23,1,-1]];
If[rep==="B2",matele=If[i===1 || i===24 || i===22 || i===23,1,-1];If[inv===1,matele=-matele]]];
(* [0,n,n] *)
(* 1:0, 16:\[Pi]; C: 22:0, 17:\[Pi] *)
If[group==="C2v",
If[rep==="A1",matele=1];
If[rep==="A2",matele=If[inv===0,1,-1]];
If[rep==="B1",matele=If[i===1 || i===22,1,-1]];
If[rep==="B2",matele=If[i===1 || i===22,1,-1];If[inv===1,matele=-matele]]];
(* [n,n,n] *)
(* 1:0, 3:(2\[Pi])/3, 2:(4\[Pi])/3; Subscript[C, xz]: 19:0, 21:(2\[Pi])/3, 17:(4\[Pi])/3 *)
If[group==="C3v",
If[rep==="A1",matele=1];
If[rep==="A2",matele=If[inv===0,1,-1]];
If[rep==="E2",If[inv===0,
If[i===1,mat={{1,0},{0,1}}];
If[i===3,mat={{-(1/2),-((I Sqrt[3])/2)},{-((I Sqrt[3])/2),-(1/2)}}];
If[i===2,mat={{-(1/2),(I Sqrt[3])/2},{(I Sqrt[3])/2,-(1/2)}}];
matele=mat[[a,b]],
If[i===19,mat={{-1,0},{0,1}}];
If[i===21,mat={{1/2,-((I Sqrt[3])/2)},{(I Sqrt[3])/2,-(1/2)}}];
If[i===17,mat={{1/2,(I Sqrt[3])/2},{-((I Sqrt[3])/2),-(1/2)}}];
matele=mat[[a,b]]]]];
(* [n,m,0] and [n,n,m] ]*)
(* 1:0; Subscript[C, z]: 24:\[Pi] and 1:0; Subscript[C, xy]: 17:\[Pi] *)
If[group==="C2nm0" || group==="C2nnm",
If[rep==="A",matele=1];
If[rep==="B",matele=If[inv===0,1,-1]]];
(* [n,m,p] ]*)
(* 1:0 *)
If[group==="C1",
If[rep==="A",matele=1]];
If[matele===10000,Throw["Wrong irrep."]];
Return[matele];
];


LittlerepD[group_,rep_,inv_,i_,a_,b_]:=Module[{X,mat,matele},
X[1]=PauliMatrix[4];
X[2]=-(1/2)PauliMatrix[4]+I Sqrt[3]/2 PauliMatrix[2];
X[3]=-(1/2)PauliMatrix[4]-I Sqrt[3]/2 PauliMatrix[2];
X[4]=-(1/2)PauliMatrix[3]-Sqrt[3]/2 PauliMatrix[1];
X[5]=PauliMatrix[3];
X[6]=-(1/2)PauliMatrix[3]+Sqrt[3]/2 PauliMatrix[1];
X[7]=I 1/Sqrt[2] (PauliMatrix[1]+PauliMatrix[2]);
X[8]=1/Sqrt[2] (PauliMatrix[1]-PauliMatrix[2]);
matele=10000;
(* [0,0,n] *)
If[group==="C4vD",
If[rep==="A1",matele=1];
If[rep==="A2",matele=If[MemberQ[{1,48,24,15,14,38,39,25},i],1,-1]];
If[rep==="B1",matele=If[MemberQ[{1,48,24,46,47,22,23,25},i],1,-1]];
If[rep==="B2",matele=If[MemberQ[{1,48,24,42,43,18,19,25},i],1,-1]];
If[rep==="E2",
If[MemberQ[{1,25},i],mat=PauliMatrix[4]];
If[MemberQ[{48,24},i],mat=-PauliMatrix[4]];
If[MemberQ[{15,39},i],mat=I PauliMatrix[3]];
If[MemberQ[{14,38},i],mat=-I PauliMatrix[3]];
If[MemberQ[{46,22},i],mat=PauliMatrix[1]];
If[MemberQ[{47,23},i],mat=-PauliMatrix[1]];
If[MemberQ[{42,18},i],mat=-PauliMatrix[2]];
If[MemberQ[{43,19},i],mat=PauliMatrix[2]];
matele=mat[[a,b]]];
If[rep==="G1" || rep==="G2",
matele=ODrep["G1",i,a,b];
If[rep==="G1",If[MemberQ[{46,47,22,23,42,43,18,19},i],matele=-matele]];
If[rep==="G2",If[MemberQ[{15,14,38,39,42,43,18,19},i],matele=-matele]]];
];
(* [0,n,n] *)
If[group==="C2vD",
If[rep==="A1",matele=1];
If[rep==="A2",matele=If[MemberQ[{1,40,16,25},i],1,-1]];
If[rep==="B1",matele=If[MemberQ[{1,46,22,25},i],1,-1]];
If[rep==="B2",matele=If[MemberQ[{1,41,17,25},i],1,-1]];
If[rep==="G1",matele=ODrep["G1",i,a,b];If[MemberQ[{41,17,46,22},i],1,-1]];
];
(* [n,n,n] *)
If[group==="C3vD",
If[rep==="A1",matele=1];
If[rep==="A2",matele=If[MemberQ[{1,3,2,26,27,25},i],1,-1]];
If[rep==="K1",
If[MemberQ[{1,26,27},i],matele=1];
If[MemberQ[{3,2,25},i],matele=-1];
If[MemberQ[{41,19,21},i],matele=I];
If[MemberQ[{43,45,17},i],matele=-I]];
If[rep==="K2",
If[MemberQ[{1,26,27},i],matele=1];
If[MemberQ[{3,2,25},i],matele=-1];
If[MemberQ[{43,45,17},i],matele=I];
If[MemberQ[{41,19,21},i],matele=-I]];
If[rep==="E2",
If[MemberQ[{1,25},i],mat=X[1]];
If[MemberQ[{3,27},i],mat=X[3]];
If[MemberQ[{2,26},i],mat=X[2]];
If[MemberQ[{41,17},i],mat=-X[4]];
If[MemberQ[{19,43},i],mat=-X[5]];
If[MemberQ[{21,45},i],mat=-X[6]];
matele=mat[[a,b]]]
];
(* [n,m,0] *)
If[group==="C2nm0D",
If[rep==="A1",matele=1];
If[rep==="A2",matele=If[MemberQ[{1,25},i],1,-1]];
If[rep==="K1",
If[i===1,matele=1];
If[i===25,matele=-1];
If[i===48,matele=I];
If[i===24,i],matele=-I];
If[rep==="K2",
If[i===1,matele=1];
If[i===25,matele=-1];
If[i===24,matele=I];
If[i===48,i],matele=-I]];
(* [n,n,m] *)
If[group==="C2nnmD",
If[rep==="A1",matele=1];
If[rep==="A2",matele=If[MemberQ[{1,25},i],1,-1]];
If[rep==="K1",
If[i===1,matele=1];
If[i===25,matele=-1];
If[i===41,matele=I];
If[i===17,i],matele=-I];
If[rep==="K2",
If[i===1,matele=1];
If[i===25,matele=-1];
If[i===17,matele=I];
If[i===41,i],matele=-I]];
(* [n,n,p] *)
If[group==="C1D",
If[rep==="A",matele=1];
If[rep==="K",matele=If[i===1,i],1,-1]];
If[matele===10000,Throw["Wrong irrep."]];
Return[matele];
];


(* Combine all irreps. *)
(* Note that for the convenience of coding, here lists only irrep. of O group instead of Subscript[O, h] . But for little groups, it lists the rep . of elements that are not pure rotations, which connects to the O group after multiplying the inversion. *)
Representation[group_,rep_,inv_,i_,a_,b_]:=Module[{},
(* Distribute the irrep. functions *)
If[MemberQ[{"Oh"},group],Return[Orep[rep,i,a,b]]];
If[MemberQ[{"C4v","C2v","C3v","C2nm0","C2nnm","C1"},group],Return[Littlerep[group,rep,inv,i,a,b]]];
If[MemberQ[{"OhD"},group],Return[ODrep[rep,i,a,b]]];
If[MemberQ[{"C4vD","C2vD","C3vD","C2nm0D","C2nnmD","C1D"},group],Return[LittlerepD[group,rep,inv,i,a,b]]];
Throw["No such group/little group."];
Return[0];
];


End[];
