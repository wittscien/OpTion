(* ::Package:: *)

subduction::usage = "subduction matrices.";
subductionLittle::usage = "subductionLittle[group,rep,r,\[Eta]tilde,\[Lambda]] gives the subduction matrices for little gropus.";


Begin["`Subduction`"];


(* Subduction *)
subduction=Association[
{0,"A1"}->{{1}},
{1,"T1"}->{{-(1/Sqrt[2]),0,1/Sqrt[2]},{1/Sqrt[2],0,1/Sqrt[2]},{0,1,0}},
{2,"T2"}->{{0,1/Sqrt[2],0,1/Sqrt[2],0},{0,1/Sqrt[2],0,-(1/Sqrt[2]),0},{1/Sqrt[2],0,0,0,-(1/Sqrt[2])}},
{2,"E"}->{{0,0,1,0,0},{1/Sqrt[2],0,0,0,1/Sqrt[2]}},
{1/2,"G1"}->{{1,0},{0,1}},
{3/2,"H"}->{{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}},
{5/2,"H"}->{{0,Sqrt[1/6],0,0,0,Sqrt[5/6]},{0,0,-1,0,0,0},{0,0,0,1,0,0},{-Sqrt[(5/6)],0,0,0,-Sqrt[(1/6)],0}},
{5/2,"G2"}->{{Sqrt[1/6],0,0,0,-Sqrt[(5/6)],0},{0,-Sqrt[(5/6)],0,0,0,Sqrt[1/6]}}
];


(* Subduction for little groups, used in one-hadron operator constructions *)
subductionLittle[group_,rep_,r_,\[Eta]tilde_,\[Lambda]_]:=Module[{S,indices},
indices={group,rep,r,Abs[\[Lambda]]};
If[MemberQ[{{"C4v","A1",1,0},{"C2v","A1",1,0},{"C3v","A1",1,0}},indices] && \[Eta]tilde===1,S=1,
If[MemberQ[{{"C4v","A2",1,0},{"C2v","A2",1,0},{"C3v","A2",1,0}},indices] && \[Eta]tilde===-1,S=1,
If[MemberQ[{{"C4v","E2",1,1},{"C4v","B1",1,2},{"C4v","E2",1,3},{"C4v","A1",1,4},{"C2v","B1",1,1},{"C2v","A1",1,2},{"C2v","B1",1,3},{"C2v","A1",1,4},{"C3v","E2",1,1},{"C3v","A2",1,3},{"C3v","E2",2,4}},indices],S=(KroneckerDelta[Sign[\[Lambda]],1]+\[Eta]tilde KroneckerDelta[Sign[\[Lambda]],-1])/Sqrt[2],
If[MemberQ[{{"C4v","E2",2,1},{"C4v","B2",1,2},{"C4v","A2",1,4},{"C2v","B2",1,1},{"C2v","A2",1,2},{"C2v","B2",1,3},{"C2v","A2",1,4},{"C3v","E2",2,1},{"C3v","E2",1,2},{"C3v","A1",1,3},{"C3v","E2",1,4}},indices],S=(KroneckerDelta[Sign[\[Lambda]],1]-\[Eta]tilde KroneckerDelta[Sign[\[Lambda]],-1])/Sqrt[2],
If[MemberQ[{{"C4v","E2",2,3}},indices],S=(-KroneckerDelta[Sign[\[Lambda]],1]+\[Eta]tilde KroneckerDelta[Sign[\[Lambda]],-1])/Sqrt[2],
If[MemberQ[{{"C3v","E2",2,2}},indices],S=(-KroneckerDelta[Sign[\[Lambda]],1]-\[Eta]tilde KroneckerDelta[Sign[\[Lambda]],-1])/Sqrt[2],
S=0;(*Print["Wrong subduction"]*)]]]]]];
Return[S];
];


End[];
