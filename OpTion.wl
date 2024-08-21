(* ::Package:: *)

(*
Author: Haobo Yan
Description: OpTion is a Mathematica package for Operator construcTion in lattice QCD.
Contact: haobo@stu.pku.edu.cn
*)


(* History *)
(* 2024.08.19: Remove linearly-dependent operators. This happens expeciaaly for N-N scatterings or high back-to-back momenta. *)
(* 2024.08.19: Add general momenta. *)
(* 2023.11.30: Stable. Build the package OpTion and release to CLQCD. *)
(* 2023.11.29: Change waverules from string to var. *)
(* 2023.11.26: Start the construction for the little group of the double cover groups. *)
(* 2023.11.25: Fix the double cover groups. The Wigner matrix was wrong *)
(* 2023.11.25: Change the ordering of the double cover groups. *)
(* 2023.11.13: Clean the code and test the double cover groups. *)
(* 2023.04.16: Add several testings to check if they transform correctly. *)
(* 2023.03.08: Fix the little groups. *)
(* 2023.03.07: Try to fix the little groups. *)
(* 2023.02.15: Add the projection method for little groups. *)
(* 2022.12.07: Fix the projection method. Everything works fine by now. *)
(* 2022.12.03: Fix the partial-wave method. Add the parity interface. *)
(* 2022.12.02: Fix the projection method. *)
(* 2022.11.30: First version. The projection method does not consistent with the reference. *)


Print[Style["OpTion: package for Operator construcTion in lattice QCD.",Bold]];
Print["Authors: Haobo Yan"];
Print["Contact: haobo@stu.pku.edu.cn"];
BeginPackage["OpTion`"];


(* Load the NC package *)
AppendTo[$Path,FileNameJoin[{DirectoryName[$InputFileName],"OpTion"}]];
<< NC`;
<< NCAlgebra`;


(* Load the package *)
(*files=Flatten[FileNames["*.wl",FileNameJoin[{DirectoryName[$InputFileName],"OpTion",#}],Infinity]&/@{"Algebra","Groups","OneHadronGeneration","TwoHadronGeneration"}];*)
files={
"NC/init.m",
"Algebra/Algebra.wl",
"Misc/Misc.wl",
"Groups/Subduction.wl",
"Groups/GroupsElements.wl",
"Groups/WignerD.wl",
"Groups/Representations.wl",
"OneHadron/Projection.wl",
"OneHadron/OneHadronOperator.wl",
"TwoHadron/Projection.wl",
"TwoHadron/PartialWave.wl",
"TwoHadron/TwoHadronOperator.wl",
"NHadron/Projection.wl",
"NHadron/NHadronOperator.wl"
};
files=Table[FileNameJoin[{DirectoryName[$InputFileName],"OpTion",file}],{file,files}];
Get/@files;
Remove[files]
EndPackage[];
