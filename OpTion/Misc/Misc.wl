(* ::Package:: *)

MomToGroup::usage = "translates the total momentum into the (little) group.";
ParToSpinParity::usage = "translate the particle labels into the spin and parity.";
Gele::usage = "the Little group list.";
parparity::usage = "parparity[parity] translate +, - to 1, -1.";
GenerateMomVectors::usage = "GenerateMomVectors[N] generate momentum vectors whose norms are less than N."
GenerateMomVectorsM::usage = "GenerateMomVectorsM[N,M,ptot] generate momentum vectors whose norms are less than N for M hadrons."
IsLinearlyIndependent::usage = "IsLinearlyIndependent[exprList,newExpr] judges if newExpr is independent on expressions in exprList."


Begin["`Misc`"];


(* Translate the total momentum into the (little) group *)
MomToGroup[ptot_,IfD_]:=Module[{x,y,z,group},
If[ptot==={0,0,0},Return["Oh"]];
Assert[Length[ptot]===3];
{x,y,z}=ptot;
group="";
If[x===y===z&&x===0,group="Oh"];
If[x===y===0&&z>0,group="C4v"];
If[x===0&&y===z&&z>0,group="C2v"];
If[x===y===z&&x>0,group="C3v"];
If[x=!=y&&x>0&&y>0&&z===0,group="C2nm0"];
If[x===y&&x>0&&z>0&&x=!=z,group="C2nnm"];
If[x>0&&y>0&&z>0&&x=!=y&&y=!=z&&x=!=z,group="C1"];
If[group==="",Throw["Change the symmetry axis."]];
If[IfD===1,group=group<>"D"];
Return[group];
];
MomToGroup[ptot_]:=MomToGroup[ptot,0];


(* Translate the particle labels into the spin and parity *)
ParToSpinParity=<|
"P"->{0,"-"},"Psudoscalar"->{0,"-"},
"S"->{0,"+"},"Scalar"->{0,"+"},
"V"->{1,"-"},"Vector"->{1,"-"},
"A"->{1,"+"},"Axialvector"->{1,"+"},
"N"->{1/2,"+"},"Nucleon"->{1/2,"+"},
"M"->{1/2,"-"},"Nucleon-"->{1/2,"-"}
|>;


(* The Little group list *)
Gele=<|
"Oh"-><|"+"->Range[24],"-"->Range[24]|>,
"C4v"-><|"+"->{1,24,14,15},"-"->{22,23,19,18}|>,
"C2v"-><|"+"->{1,16},"-"->{22,17}|>,
"C3v"-><|"+"->{1,3,2},"-"->{19,21,17}|>,
"C2nm0"-><|"+"->{1},"-"->{24}|>,
"C2nnm"-><|"+"->{1},"-"->{19}|>,
"C1"-><|"+"->{1},"-"->{}|>,
"OhD"-><|"+"->Range[48],"-"->Range[48]|>,
"C4vD"-><|"+"->{1,48,24,15,14,38,39,25},"-"->{46,47,22,23,42,43,18,19}|>,
"C2vD"-><|"+"->{1,40,16,25},"-"->{41,17,46,22}|>,
"C3vD"-><|"+"->{1,3,2,26,27,25},"-"->{41,19,21,43,45,17}|>,
"C2nm0D"-><|"+"->{1,25},"-"->{48,24}|>,
"C2nnmD"-><|"+"->{1,25},"-"->{43,19}|>,
"C1D"-><|"+"->{1,25},"-"->{}|>
|>;


parparity[parity_]:=If[parity==="+",Return[1],If[parity==="-",Return[-1],Print["wrong particle parity"]]];


(* Generate momentum vectors whose norms are less than N *)
GenerateMomVectors[N_]:=Module[{},
If[!IntegerQ[N] || N<0,Throw["N should be an positive integer.","f"]];Return[SortBy[Select[Tuples[Range[-N,N],3],Norm[#]<=N&],Norm]]];


(* For M-hadron operators *)
GenerateMomVectorsM[N_,M_,ptot_]:=Module[{},
Return[SortBy[Select[Tuples[GenerateMomVectors[N],M],Total[#]===ptot &],Total[(Norm/@#)^2]&]]];


(* 2024.08.19: Check linear dependency *)
IsLinearlyIndependent[exprList_,newExpr_]:=Module[{allMonomials,coeffMatrix,augmentedMatrix},
(* MatrixRank refuses the empty matrix *)
If[exprList==={},Return[True]];
(* Extract all distinct monomials from exprList and newExpr *)
allMonomials=Union[ONormalize/@Flatten[Map[List@@Expand[#]&,Join[exprList,{newExpr}]]/. Plus->List]];
(* Compute the coefficient matrix for the list with respect to the monomials *)
coeffMatrix=Table[Coefficient[expr,mon],{expr,exprList},{mon,allMonomials}];
(* Compute the coefficients for the new expression with respect to the monomials *)
augmentedMatrix=Append[coeffMatrix,Coefficient[newExpr,#]&/@allMonomials];
(* Check if the new expression introduces a linear dependence *)
Return[MatrixRank[augmentedMatrix]>MatrixRank[coeffMatrix]];
]


End[];
