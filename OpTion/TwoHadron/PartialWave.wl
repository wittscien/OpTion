(* ::Package:: *)

PartialWaveOperatorProjected::usage = "PartialWaveOperatorProjected[rep,r,p,J,L,S,s1,p1,s2,p2] generates two-hadron operators by the partial wave method.";


Begin["`Projection`"];


(* Partial-wave method *)
PartialWaveOperator[p_,J_,mJ_,L_,S_,s1_,p1_,s2_,p2_]:=Module[{newp,ope,\[Theta],\[Phi],beforeprojection},
ope=0;
(* 2025.06.02: This is equivalent to define Y_l^m = 0 for |p|!=0 and L!=0. The line "If[p==={0,0,0} && L!=0,Return[0]];" was commented out but I forgot why. Now I recover this line. *)
(* 2023.11.30: The comment last year seems to work only for mesons. But I am curious how baryon, without relative momenta, sustains orbital angular momentum. *)
(* 2022.12.03: If |p|==0, there is no way to have angular momentum. So PV |p|=0 Subscript[A, 1] rep. has no operator. The Eq. (4.11) in Prelovsek's paper should be refined. *)
If[p==={0,0,0} && L!=0,Return[0]];
Do[Do[Do[Do[If[mL+mS===mJ && ms1+ms2===mS,
beforeprojection=0;
Do[newp=RotationMatrix[Oh["\[Omega]"][[i]],Oh["n"][[i]]] . p;
If[newp[[1]]!=0 || newp[[2]]!=0,
(* \[Theta]=ArcCos[newp[[3]]/Norm[newp]];
\[Phi]=Sign[newp[[2]]]ArcCos[newp[[1]]/Sqrt[newp[[1]]^2+newp[[2]]^2]]; *)
(* Should be very careful at special cases when only z!=0 *)
\[Theta]=ToSphericalCoordinates[newp][[2]];
\[Phi]=ToSphericalCoordinates[newp][[3]];,If[newp[[3]]>0,\[Theta]=0;\[Phi]=0,\[Theta]=\[Pi];\[Phi]=0]];
beforeprojection+=Conjugate[SphericalHarmonicY[L,mL,\[Theta],\[Phi]]]H[1,s1,ms1,p1][newp]H[2,s2,ms2,p2][-newp],{i,1,24}];
ope+=ClebschGordan[{L,mL},{S,mS},{J,mJ}] ClebschGordan[{s1,ms1},{s2,ms2},{S,mS}] beforeprojection],{ms1,-s1,s1}],{ms2,-s2,s2}],{mS,-S,S}],{mL,-L,L}];
Return[ope];
];


(* Partial-wave method *)
PartialWaveOperatorD[p_,J_,mJ_,L_,S_,s1_,p1_,s2_,p2_]:=Module[{newp,ope,\[Theta],\[Phi],beforeprojection},
ope=0;
If[p==={0,0,0} && L!=0,Return[0]];
Do[Do[Do[Do[If[mL+mS===mJ && ms1+ms2===mS,
beforeprojection=0;
Do[newp=RotationMatrix[OhD["\[Omega]"][[i]],OhD["n"][[i]]] . p;
If[newp[[1]]!=0 || newp[[2]]!=0,
(* \[Theta]=ArcCos[newp[[3]]/Norm[newp]];
\[Phi]=Sign[newp[[2]]]ArcCos[newp[[1]]/Sqrt[newp[[1]]^2+newp[[2]]^2]]; *)
(* Should be very careful at special cases when only z!=0 *)
\[Theta]=ToSphericalCoordinates[newp][[2]];
\[Phi]=ToSphericalCoordinates[newp][[3]];,If[newp[[3]]>0,\[Theta]=0;\[Phi]=0,\[Theta]=\[Pi];\[Phi]=0]];
beforeprojection+=Conjugate[SphericalHarmonicY[L,mL,\[Theta],\[Phi]]]H[1,s1,ms1,p1][newp]H[2,s2,ms2,p2][-newp],{i,1,48}];
ope+=ClebschGordan[{L,mL},{S,mS},{J,mJ}] ClebschGordan[{s1,ms1},{s2,ms2},{S,mS}] beforeprojection],{ms1,-s1,s1}],{ms2,-s2,s2}],{mS,-S,S}],{mL,-L,L}];
Return[ope];
];


PartialWaveOperatorProjected[rep_,r_,p_,J_,L_,S_,s1_,p1_,s2_,p2_]:=Module[{ope},
ope=0;
If[(rep==="A1" && MemberQ[{0},J]) || (rep==="T1" && MemberQ[{1,3,4},J]) || (rep==="T2" && MemberQ[{2,3,4},J]) || (rep==="E" && MemberQ[{2,4},J]) || (rep==="G1" && MemberQ[{1/2},J]) || (rep==="H" && MemberQ[{3/2,5/2},J]) || (rep==="G2" && MemberQ[{5/2},J]),
Do[ope+=subduction[{J,rep}][[r,-mJ+J+1]]PartialWaveOperator[p,J,mJ,L,S,s1,p1,s2,p2],{mJ,-J,J}],
If[(rep==="G1" && MemberQ[{1/2},J]) || (rep==="H" && MemberQ[{3/2,5/2},J]) || (rep==="G2" && MemberQ[{5/2},J]),
Do[ope+=subduction[{J,rep}][[r,-mJ+J+1]]PartialWaveOperatorD[p,J,mJ,L,S,s1,p1,s2,p2],{mJ,-J,J}]]];
ope=OSimplify[ope];
Return[ope];
];


End[];
