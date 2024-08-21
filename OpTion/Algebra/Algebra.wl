(* ::Package:: *)

symbols::usage = "replacing lists for momentum";
(*waverules::usage = "replacing lists for two-hadron string operators";*)
waverulesvar::usage = "replacing lists for two-hadron variable operators";
waverulesOne::usage = "replacing lists for one-hadron operators";
varToStringOne::usage = "replacing lists for one-hadron operators";
PauliMatrix3d2::usage = "PauliMatrix3d2[k] gives the \!\(\*SuperscriptBox[\(k\), \(th\)]\) generalized Pauli matrix \!\(\*SubscriptBox[\(\[Sigma]\), \(k\)]\) for spin s=\!\(\*FractionBox[\(3\), \(2\)]\).";
ONormalize::usage = "ONormalize[ex] discard the coefficients of an expression. Also used in my contraction tools.";
OSimplify::usage = "OSimplify[ex] transforms the operator expression from its spherical forms into Cartesian vector forms and normalizes the expression.";
NCDSimplify::usage = "NCDSimplify[ex] transforms the operator expression from its spherical forms into Cartesian vector forms and normalizes the expression.";


Begin["`Algebra`"];


(* Just notations *)
symbols={
{0,0,0}->0,
{1,0,0}->Subscript["e", "x"],{-1,0,0}->-Subscript["e", "x"],{0,1,0}->Subscript["e", "y"],{0,-1,0}->-Subscript["e", "y"],{0,0,1}->Subscript["e", "z"],{0,0,-1}->-Subscript["e", "z"],
{0,1,1}->Subscript["e", "y","z"],{1,0,1}->Subscript["e", "x","z"],{1,1,0}->Subscript["e", "x","y"],{0,-1,1}->Subscript["e", -"y","z"],{-1,0,1}->Subscript["e", -"x","z"],{-1,1,0}->Subscript["e", -"x","y"],{0,1,-1}->Subscript["e", "y",-"z"],{1,0,-1}->Subscript["e", "x",-"z"],{1,-1,0}->Subscript["e", "x",-"y"],{0,-1,-1}->Subscript["e", -"y",-"z"],{-1,0,-1}->Subscript["e", -"x",-"z"],{-1,-1,0}->Subscript["e", -"x",-"y"],
{1,1,1}->Subscript["e", "x","y","z"],{-1,1,1}->Subscript["e", -"x","y","z"],{1,-1,1}->Subscript["e", "x",-"y","z"],{1,1,-1}->Subscript["e", "x","y",-"z"],{-1,-1,1}->Subscript["e", -"x",-"y","z"],{-1,1,-1}->Subscript["e", -"x","y",-"z"],{1,-1,-1}->Subscript["e", "x",-"y",-"z"],{-1,-1,-1}->Subscript["e", -"x",-"y",-"z"],
{0,0,2}->Subscript["e", 2"z"],{0,2,0}->Subscript["e", 2"y"],{2,0,0}->Subscript["e", 2"x"],{0,0,-2}->Subscript["e", -2"z"],{0,-2,0}->Subscript["e", -2"y"],{-2,0,0}->Subscript["e", -2"x"]
};

waverules={
\!\(\*SubscriptBox[\("\<H1\>"\), \({0, "\<-\>"}\)]\)[0][p_]->"\!\(\*SubscriptBox[\(P\), \(1\)]\)"[p],\!\(\*SubscriptBox[\("\<H2\>"\), \({0, "\<-\>"}\)]\)[0][p_]->"\!\(\*SubscriptBox[\(P\), \(2\)]\)"[p],
\!\(\*SubscriptBox[\("\<H1\>"\), \({0, "\<+\>"}\)]\)[0][p_]->"\!\(\*SubscriptBox[\(S\), \(1\)]\)"[p],\!\(\*SubscriptBox[\("\<H2\>"\), \({0, "\<+\>"}\)]\)[0][p_]->"\!\(\*SubscriptBox[\(S\), \(2\)]\)"[p],
\!\(\*SubscriptBox[\("\<H1\>"\), \({1, "\<-\>"}\)]\)[1][p_]->(-Subscript["\!\(\*SubscriptBox[\(V\), \(1\)]\)", "x"][p]+I Subscript["\!\(\*SubscriptBox[\(V\), \(1\)]\)", "y"][p])/Sqrt[2],\!\(\*SubscriptBox[\("\<H1\>"\), \({1, "\<-\>"}\)]\)[-1][p_]->(Subscript["\!\(\*SubscriptBox[\(V\), \(1\)]\)", "x"][p]+I Subscript["\!\(\*SubscriptBox[\(V\), \(1\)]\)", "y"][p])/Sqrt[2],\!\(\*SubscriptBox[\("\<H1\>"\), \({1, "\<-\>"}\)]\)[0][p_]->Subscript["\!\(\*SubscriptBox[\(V\), \(1\)]\)", "z"][p],
\!\(\*SubscriptBox[\("\<H2\>"\), \({1, "\<-\>"}\)]\)[1][p_]->(-Subscript["\!\(\*SubscriptBox[\(V\), \(2\)]\)", "x"][p]+I Subscript["\!\(\*SubscriptBox[\(V\), \(2\)]\)", "y"][p])/Sqrt[2],\!\(\*SubscriptBox[\("\<H2\>"\), \({1, "\<-\>"}\)]\)[-1][p_]->(Subscript["\!\(\*SubscriptBox[\(V\), \(2\)]\)", "x"][p]+I Subscript["\!\(\*SubscriptBox[\(V\), \(2\)]\)", "y"][p])/Sqrt[2],\!\(\*SubscriptBox[\("\<H2\>"\), \({1, "\<-\>"}\)]\)[0][p_]->Subscript["\!\(\*SubscriptBox[\(V\), \(2\)]\)", "z"][p],
\!\(\*SubscriptBox[\("\<H1\>"\), \({1, "\<+\>"}\)]\)[1][p_]->(-Subscript["\!\(\*SubscriptBox[\(A\), \(1\)]\)", "x"][p]+I Subscript["\!\(\*SubscriptBox[\(A\), \(1\)]\)", "y"][p])/Sqrt[2],\!\(\*SubscriptBox[\("\<H1\>"\), \({1, "\<+\>"}\)]\)[-1][p_]->(Subscript["\!\(\*SubscriptBox[\(A\), \(1\)]\)", "x"][p]+I Subscript["\!\(\*SubscriptBox[\(A\), \(1\)]\)", "y"][p])/Sqrt[2],\!\(\*SubscriptBox[\("\<H1\>"\), \({1, "\<+\>"}\)]\)[0][p_]->Subscript["\!\(\*SubscriptBox[\(A\), \(1\)]\)", "z"][p],
\!\(\*SubscriptBox[\("\<H2\>"\), \({1, "\<+\>"}\)]\)[1][p_]->(-Subscript["\!\(\*SubscriptBox[\(A\), \(2\)]\)", "x"][p]+I Subscript["\!\(\*SubscriptBox[\(A\), \(2\)]\)", "y"][p])/Sqrt[2],\!\(\*SubscriptBox[\("\<H2\>"\), \({1, "\<+\>"}\)]\)[-1][p_]->(Subscript["\!\(\*SubscriptBox[\(A\), \(2\)]\)", "x"][p]+I Subscript["\!\(\*SubscriptBox[\(A\), \(2\)]\)", "y"][p])/Sqrt[2],\!\(\*SubscriptBox[\("\<H2\>"\), \({1, "\<+\>"}\)]\)[0][p_]->Subscript["\!\(\*SubscriptBox[\(A\), \(2\)]\)", "z"][p],
\!\(\*SubscriptBox[\("\<H1\>"\), \({
\*FractionBox[\(1\), \(2\)], "\<+\>"}\)]\)[1/2][p_]->Subscript["\!\(\*SubscriptBox[\(N\), \(1\)]\)", (1/2)][p],\!\(\*SubscriptBox[\("\<H1\>"\), \({
\*FractionBox[\(1\), \(2\)], "\<+\>"}\)]\)[-(1/2)][p_]->Subscript["\!\(\*SubscriptBox[\(N\), \(1\)]\)", -(1/2)][p],
\!\(\*SubscriptBox[\("\<H2\>"\), \({
\*FractionBox[\(1\), \(2\)], "\<+\>"}\)]\)[1/2][p_]->Subscript["\!\(\*SubscriptBox[\(N\), \(2\)]\)", (1/2)][p],\!\(\*SubscriptBox[\("\<H2\>"\), \({
\*FractionBox[\(1\), \(2\)], "\<+\>"}\)]\)[-(1/2)][p_]->Subscript["\!\(\*SubscriptBox[\(N\), \(2\)]\)", -(1/2)][p],
\!\(\*SubscriptBox[\("\<H1\>"\), \({
\*FractionBox[\(1\), \(2\)], "\<-\>"}\)]\)[1/2][p_]->Subscript["\!\(\*SubscriptBox[\(M\), \(1\)]\)", (1/2)][p],\!\(\*SubscriptBox[\("\<H1\>"\), \({
\*FractionBox[\(1\), \(2\)], "\<-\>"}\)]\)[-(1/2)][p_]->Subscript["\!\(\*SubscriptBox[\(M\), \(1\)]\)", -(1/2)][p],
\!\(\*SubscriptBox[\("\<H2\>"\), \({
\*FractionBox[\(1\), \(2\)], "\<-\>"}\)]\)[1/2][p_]->Subscript["\!\(\*SubscriptBox[\(M\), \(2\)]\)", (1/2)][p],\!\(\*SubscriptBox[\("\<H2\>"\), \({
\*FractionBox[\(1\), \(2\)], "\<-\>"}\)]\)[-(1/2)][p_]->Subscript["\!\(\*SubscriptBox[\(M\), \(2\)]\)", -(1/2)][p]
};

waverulesvar={
OpTion`Projection`H[i_,0,0,"-"][p_]:>("\!\(\*SubscriptBox[\(P\), \("<>ToString[i]<>"\)]\)")[p],
OpTion`Projection`H[i_,0,0,"+"][p_]:>("\!\(\*SubscriptBox[\(S\), \("<>ToString[i]<>"\)]\)")[p],
OpTion`Projection`H[i_,1,1,"-"][p_]:>(-Subscript["\!\(\*SubscriptBox[\(V\), \("<>ToString[i]<>"\)]\)", "x"][p]+I Subscript["\!\(\*SubscriptBox[\(V\), \("<>ToString[i]<>"\)]\)", "y"][p])/Sqrt[2],
OpTion`Projection`H[i_,1,-1,"-"][p_]:>(Subscript["\!\(\*SubscriptBox[\(V\), \("<>ToString[i]<>"\)]\)", "x"][p]+I Subscript["\!\(\*SubscriptBox[\(V\), \("<>ToString[i]<>"\)]\)", "y"][p])/Sqrt[2],
OpTion`Projection`H[i_,1,0,"-"][p_]:>Subscript["\!\(\*SubscriptBox[\(V\), \("<>ToString[i]<>"\)]\)", "z"][p],
OpTion`Projection`H[i_,1,1,"+"][p_]:>(-Subscript["\!\(\*SubscriptBox[\(A\), \("<>ToString[i]<>"\)]\)", "x"][p]+I Subscript["\!\(\*SubscriptBox[\(A\), \("<>ToString[i]<>"\)]\)", "y"][p])/Sqrt[2],
OpTion`Projection`H[i_,1,-1,"+"][p_]:>(Subscript["\!\(\*SubscriptBox[\(A\), \("<>ToString[i]<>"\)]\)", "x"][p]+I Subscript["\!\(\*SubscriptBox[\(A\), \("<>ToString[i]<>"\)]\)", "y"][p])/Sqrt[2],
OpTion`Projection`H[i_,1,0,"+"][p_]:>Subscript["\!\(\*SubscriptBox[\(A\), \("<>ToString[i]<>"\)]\)", "z"][p],
OpTion`Projection`H[i_,1/2,1/2,"+"][p_]:>Subscript["\!\(\*SubscriptBox[\(N\), \("<>ToString[i]<>"\)]\)", (1/2)][p],
OpTion`Projection`H[i_,1/2,-(1/2),"+"][p_]:>Subscript["\!\(\*SubscriptBox[\(N\), \("<>ToString[i]<>"\)]\)", -(1/2)][p],
OpTion`Projection`H[i_,1/2,1/2,"-"][p_]:>Subscript["\!\(\*SubscriptBox[\(M\), \("<>ToString[i]<>"\)]\)", (1/2)][p],
OpTion`Projection`H[i_,1/2,-(1/2),"-"][p_]:>Subscript["\!\(\*SubscriptBox[\(M\), \("<>ToString[i]<>"\)]\)", -(1/2)][p]
};


(* Definitions *)
waverulesOne={
Subscript[OpTion`Projection`V, 1]->-(I/Sqrt[2])(Subscript[OpTion`Projection`V, "x"]-I Subscript[OpTion`Projection`V, "y"]),Subscript[OpTion`Projection`V, -1]->I/Sqrt[2] (Subscript[OpTion`Projection`V, "x"]+I Subscript[OpTion`Projection`V, "y"]),Subscript[OpTion`Projection`V, 0]->I Subscript[OpTion`Projection`V, "z"],
Subscript[OpTion`Projection`A, 1]->-(I/Sqrt[2])(Subscript[OpTion`Projection`A, "x"]-I Subscript[OpTion`Projection`A, "y"]),Subscript[OpTion`Projection`A, -1]->I/Sqrt[2] (Subscript[OpTion`Projection`A, "x"]+I Subscript[OpTion`Projection`A, "y"]),Subscript[OpTion`Projection`A, 0]->I Subscript[OpTion`Projection`A, "z"],
Subscript[D, 1]->-(I/Sqrt[2])(Subscript[D, "x"]-I Subscript[D, "y"]),Subscript[D, -1]->I/Sqrt[2] (Subscript[D, "x"]+I Subscript[D, "y"]),Subscript[D, 0]->I Subscript[D, "z"]};

varToStringOne={OpTion`Projection`P->"P",OpTion`Projection`S->"S",OpTion`Projection`V->"V",OpTion`Projection`A->"A"};


(* Pauli matrix for spin-3/2 *)
PauliMatrix3d2[1]=\!\(\*
TagBox[
RowBox[{"(", GridBox[{
{"0", 
SqrtBox["3"], "0", "0"},
{
SqrtBox["3"], "0", "2", "0"},
{"0", "2", "0", 
SqrtBox["3"]},
{"0", "0", 
SqrtBox["3"], "0"}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\);
PauliMatrix3d2[2]=I \!\(\*
TagBox[
RowBox[{"(", GridBox[{
{"0", 
RowBox[{"-", 
SqrtBox["3"]}], "0", "0"},
{
SqrtBox["3"], "0", 
RowBox[{"-", "2"}], "0"},
{"0", "2", "0", 
RowBox[{"-", 
SqrtBox["3"]}]},
{"0", "0", 
SqrtBox["3"], "0"}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.7]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\);
PauliMatrix3d2[3]=({
 {3, 0, 0, 0},
 {0, 1, 0, 0},
 {0, 0, -1, 0},
 {0, 0, 0, -3}
});
PauliMatrix3d2[4]=({
 {1, 0, 0, 0},
 {0, 1, 0, 0},
 {0, 0, 1, 0},
 {0, 0, 0, 1}
});


ONormalize[ex_]:=Module[{constants,expReturn},
If[Head[ex]===Times,
constants=Cases[ex,_?NumericQ,1];
If[Simplify[Times@@constants]=!=0,
expReturn=ex/(Times@@constants),expReturn=0],
ex]];


(* Simplify function *)
OSimplify[ex_]:=Module[{expfactor,expreim,expre,expim,expreturn},
expfactor=ex//.symbols//.waverulesvar//Factor;
expfactor=ONormalize[expfactor];
expreim=ComplexExpand[ReIm[expfactor]];
expre=expreim[[1]]//Expand;
expim=expreim[[2]]//Expand;
expreturn=expre+I expim;
expreturn=Expand[ONormalize[Factor[expreturn]]];
Return[expreturn];
];


(* Simplify function *)
DSimplify[ex_]:=Module[{expexpand},
(* expexpand=ex//.waverulesOne//Simplify; *)
expexpand=ex//.waverulesOne//Simplify;
(* Print["----------------------------------------------------------------------------------------------------------"]; *)
Return[expexpand];
];

NCDSimplify[ex_]:=Module[{expfactor},
expfactor=ex//.waverulesOne//NCSimplifyRational//Simplify;
expfactor=ONormalize[expfactor];
expfactor=expfactor//Factor;
expfactor=expfactor//.varToStringOne;
Return[expfactor];
];


End[];
