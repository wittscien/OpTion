(* ::Package:: *)

MyD::usage = "MyD[j,m,mp,i] gives my version of the Wigner D matrices.";
MyDD::usage = "MyDD[j,m,mp,i] gives my version of the Wigner D matrices.";


Begin["`WignerD`"];


(* Wigner D matrix *)
MyD[j_,m_,mp_,i_]:=Module[{myd,F,\[Alpha],\[Beta],\[Gamma]},
myd=0;
F=If[!IntegerQ[j] && !MemberQ[{1,2,3,4,8,9,10,14,15,16,20,21,22,23,28,29,30,31,36,37,38,39,40,41},i],-1,1];
{\[Alpha],\[Beta],\[Gamma]}=EulerAngles[RotationMatrix[Oh["\[Omega]"][[i]],Oh["n"][[i]]]]//Simplify;
(* myd=F WignerD[{j,m,mp},-Oh["\[Alpha]"][[i]],-Oh["\[Beta]"][[i]],-Oh["\[Gamma]"][[i]]]; *)
myd=F WignerD[{j,m,mp},-\[Alpha],-\[Beta],-\[Gamma]];
Return[myd];
];

MyDD[j_,m_,mp_,i_]:=Module[{myd,F,\[Alpha],\[Beta],\[Gamma]},
myd=0;
F=If[!IntegerQ[j] && MemberQ[Join[{2,4,12,16,18,20,21,23,24}, {2,4,12,16,18,20,21,23,24}+24],i],-1,1];
If[!IntegerQ[j] && 25<=i<=48,F=-F];
{\[Alpha],\[Beta],\[Gamma]}=EulerAngles[RotationMatrix[OhD["\[Omega]"][[i]],OhD["n"][[i]]]]//Simplify;
myd=F WignerD[{j,m,mp},-\[Alpha],-\[Beta],-\[Gamma]];
(* The following is equavilent to the above for j=1/2 *)
(*
mm=If[m===1/2,1,2]; mmp=If[mp===1/2,1,2];
myd2=FullSimplify[MatrixExp[-(I/2)(OhD["n"][[i]][[1]] PauliMatrix[1]+OhD["n"][[i]][[2]] PauliMatrix[2]+OhD["n"][[i]][[3]] PauliMatrix[3])OhD["\[Omega]"][[i]]]][[mm,mmp]];
*)
Return[myd];
];


End[];
