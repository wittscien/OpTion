(* ::Package:: *)

subduction::usage = "subduction matrices.";
subductionLittle::usage = "subductionLittle[group,rep,r,\[Eta]tilde,\[Lambda]] gives the subduction matrices for little gropus.";


Begin["`Subduction`"];


(* Subduction *)
(* 2025.06.02: Discussion with Sasa: T1 and T2 has different basis and thus the subduction matrices from the HSC should be transformed. *)
(* [HSC] PRD 2010 Dudek; [Prelovsek] JHEP 2017 B.4 for Subscript[T, 1] and Subscript[T, 2] *)

(* Old, not good for T1 and T2 with J>=3. *)
(*subduction=Association[
{0,"A1"}->{{1}},
{1,"T1"}->{{-(1/Sqrt[2]),0,1/Sqrt[2]},{1/Sqrt[2],0,1/Sqrt[2]},{0,1,0}},
{2,"T2"}->{{0,1/Sqrt[2],0,1/Sqrt[2],0},{0,1/Sqrt[2],0,-(1/Sqrt[2]),0},{1/Sqrt[2],0,0,0,-(1/Sqrt[2])}},
{2,"E"}->{{0,0,1,0,0},{1/Sqrt[2],0,0,0,1/Sqrt[2]}},
{3,"T1"}->{{0,0,Sqrt[(3/8)],0,0,0,Sqrt[(5/8)]},{0,0,0,-1,0,0,0},{Sqrt[(5/8)],0,0,0,Sqrt[(3/8)],0,0}},
{3,"T2"}->{{0,0,Sqrt[(5/8)],0,0,0,-Sqrt[(3/8)]},{0,-(1/Sqrt[2]),0,0,0,-(1/Sqrt[2]),0},{Sqrt[(3/8)],0,0,0,-Sqrt[(5/8)],0,0}},
{3,"A2"}->{{0,1/Sqrt[2],0,0,0,-(1/Sqrt[2]),0}},
{4,"T1"}->{{0,0,0,-Sqrt[(7/8)],0,0,0,-(1/Sqrt[8]),0},{1/Sqrt[2],0,0,0,0,0,0,0,-(1/Sqrt[2])},{0,Sqrt[(1/8)],0,0,0,Sqrt[(7/8)],0,0,0}},
{4,"T2"}->{{0,0,0,-Sqrt[(1/8)],0,0,0,Sqrt[(7/8)],0},{0,0,1/Sqrt[2],0,0,0,-(1/Sqrt[2]),0,0},{0,Sqrt[(7/8)],0,0,0,-Sqrt[(1/8)],0,0,0}},
{4,"E"}->{{Sqrt[(7/24)],0,0,0,-Sqrt[(5/12)],0,0,0,Sqrt[(7/24)]},{0,0,1/Sqrt[2],0,0,0,1/Sqrt[2],0,0}},
{1/2,"G1"}->{{1,0},{0,1}},
{3/2,"H"}->{{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}},
{5/2,"H"}->{{0,Sqrt[1/6],0,0,0,Sqrt[5/6]},{0,0,-1,0,0,0},{0,0,0,1,0,0},{-Sqrt[(5/6)],0,0,0,-Sqrt[(1/6)],0}},
{5/2,"G2"}->{{Sqrt[1/6],0,0,0,-Sqrt[(5/6)],0},{0,-Sqrt[(5/6)],0,0,0,Sqrt[1/6]}}
];
*)

subduction=Association[
{0,"A1"}->{{1}},
{1,"T1"}->{{-(1/Sqrt[2]),0,1/Sqrt[2]},{1/Sqrt[2],0,1/Sqrt[2]},{0,1,0}},
{2,"T2"}->{{0,1/Sqrt[2],0,1/Sqrt[2],0},{0,1/Sqrt[2],0,-(1/Sqrt[2]),0},{1/Sqrt[2],0,0,0,-(1/Sqrt[2])}},
{2,"E"}->{{0,0,1,0,0},{1/Sqrt[2],0,0,0,1/Sqrt[2]}},
{3,"T1"}->{{Sqrt[(5/16)],0,-Sqrt[(3/16)],0,Sqrt[(3/16)],0,-Sqrt[(5/16)]},{Sqrt[(5/16)],0,Sqrt[(3/16)],0,Sqrt[(3/16)],0,Sqrt[(5/16)]},{0,0,0,-1,0,0,0}},
{3,"T2"}->{{Sqrt[(3/16)],0,Sqrt[(5/16)],0,-Sqrt[(5/16)],0,-Sqrt[(3/16)]},{-Sqrt[(3/16)],0,Sqrt[(5/16)],0,Sqrt[(5/16)],0,-Sqrt[(3/16)]},{0,-(1/Sqrt[2]),0,0,0,-(1/Sqrt[2]),0}},
{3,"A2"}->{{0,1/Sqrt[2],0,0,0,-(1/Sqrt[2]),0}},
{4,"T1"}->{{0,Sqrt[(1/16)],0,Sqrt[(7/16)],0,Sqrt[(7/16)],0,1/Sqrt[16],0},{0,Sqrt[(1/16)],0,-Sqrt[(7/16)],0,Sqrt[(7/16)],0,-(1/Sqrt[16]),0},{1/Sqrt[2],0,0,0,0,0,0,0,-(1/Sqrt[2])}},
{4,"T2"}->{{0,Sqrt[(7/16)],0,-Sqrt[(1/16)],0,-Sqrt[(1/16)],0,Sqrt[(7/16)],0},{0,-Sqrt[(7/16)],0,-Sqrt[(1/16)],0,Sqrt[(1/16)],0,Sqrt[(7/16)],0},{0,0,1/Sqrt[2],0,0,0,-(1/Sqrt[2]),0,0}},
{4,"E"}->{{Sqrt[(7/24)],0,0,0,-Sqrt[(5/12)],0,0,0,Sqrt[(7/24)]},{0,0,1/Sqrt[2],0,0,0,1/Sqrt[2],0,0}},
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
