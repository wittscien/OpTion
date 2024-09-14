(* ::Package:: *)

OneHadronOperatorAll::usage = "OneHadronOperatorAll[ptot,rep,r,MaxND] generates all possible one-hadron operators.";
OneHadronOperatorHelicity::usage = "OneHadronOperatorHelicity[ptot,rep,r,J,\[Lambda],JD,J13,TG,ND] generates one-hadron operators.";
OneHadronOperatorHelicityAll::usage = "OneHadronOperatorHelicityAll[ptot,rep,r,MaxND] generates all possible one-hadron operators using the helicity method.";


Begin["`OneHadronOperator`"];


(* 2024.9.14: Happy Moon Festival! Here I implement my way of building one-hadron operators *)
OneHadronOperatorAll[ptot_,rep_,r_,MaxND_]:=Module[{group,JG,JD,J13,NJD,repparity,repO,P,PGD,opList,op},
group=MomToGroup[ptot];
opList={};
Do[
(* JD is useful for ND>=2 *)
NJD=If[ND>=2,2,0];
(* ND>=3 is not useful for now *)
J13=0;
Do[
JG=If[TG==="S" || TG==="P",0,If[TG==="V" || TG==="A",1,-100]];
Do[Do[Do[
op=ProjectionOperator1[ptot,rep,r,J,M,JD,J13,TG,ND];
(* I commented out IsLinearlyIndependent[opList,op] since for one-hadron operators V_z will be seen as V and z in that function. *)
If[!MemberQ[opList,op] && op =!= 0,AppendTo[opList,op]],{M,-J,J}],{JD,0,NJD}],{J,Abs[JG-ND],JG+ND}],{TG,{"S","P","V","A"}}],{ND,0,MaxND}];
(* Abs[JG-ND] could change to 0 if there are two derivatives *)
Return[opList];
];


OneHadronOperatorHelicity[ptot_,rep_,r_,J_,\[Lambda]_,JD_,J13_,TG_,ND_]:=Module[{ope,group},
group=MomToGroup[ptot];
If[group==="Oh",
ope=OSub[rep,r,J,JD,J13,TG,ND],
ope=OSub[group,rep,r,J,\[Lambda],JD,J13,TG,ND]
];
Return[ope];
];


(* Moving frame: Search all possible nonzero operators *)
(* Only for MaxND <= 1 for the current verison. 2024.9.14: I don't plan to develop further for the helicity method since I have a better way. *)
OneHadronOperatorHelicityAll[ptot_,rep_,r_,MaxND_]:=Module[{group,JG,JD,J13,repparity,repO,P,PGD,opList,op},
(*Print[group,": ",rep,": ","##################################################################################################"];*)
group=MomToGroup[ptot];
opList={};
If[group==="Oh",
(* Rest frame *)
Do[Do[
JG=If[TG==="S" || TG==="P",0,If[TG==="V" || TG==="A",1,-100]];
Do[
JD=0;
J13=0;
repO=StringTake[rep,{1,-2}];
repparity=If[StringPart[rep,-1]==="+",1,If[StringPart[rep,-1]==="-",-1,Print["wrong rep. parity"]]];
P=If[TG==="S" || TG==="A",1,If[TG==="P" || TG==="V",-1,0]];
PGD=P (-1)^ND;
If[repparity===PGD,
If[(repO==="A1" && J===0) || (repO==="T1" && J===1) || (repO==="T2" && J===2) || (repO==="E" && J===2),
op=OSub[rep,r,J,JD,J13,TG,ND];
If[!MemberQ[opList,op] && op =!= 0,AppendTo[opList,op]]]],{J,Abs[JG-ND],JG+ND}],{TG,{"S","P","V","A"}}],{ND,0,MaxND}],
(* Moving frame *)
Do[Do[
JG=If[TG==="S" || TG==="P",0,If[TG==="V" || TG==="A",1,-100]];
Do[Do[
JD=0;
J13=0;
op=OSub[group,rep,r,J, \[Lambda],JD,J13,TG,ND];
If[!MemberQ[opList,op] && op =!= 0,AppendTo[opList,op]],{\[Lambda],0,J}],{J,Abs[JG-ND],JG+ND}],{TG,{"S","P","V","A"}}],{ND,0,MaxND}]];
(* Abs[JG-ND] could change to 0 if there are two derivatives *)
Return[opList];
];


OneHadronOperatorAllOld[ptot_,rep_,r_,MaxND_]:=Module[{group,JG,JD,J13,repparity,repO,P,PGD,opList,op},
(*Print[group,": ",rep,": ","##################################################################################################"];*)
group=MomToGroup[ptot];
opList={};
If[group==="Oh",
(* Rest frame *)
Do[Do[
JG=If[TG==="S" || TG==="P",0,If[TG==="V" || TG==="A",1,-100]];
Do[
JD=0;
J13=0;
repO=StringTake[rep,{1,-2}];
repparity=If[StringPart[rep,-1]==="+",1,If[StringPart[rep,-1]==="-",-1,Print["wrong rep. parity"]]];
P=If[TG==="S" || TG==="A",1,If[TG==="P" || TG==="V",-1,0]];
PGD=P (-1)^ND;
If[repparity===PGD,
If[(repO==="A1" && J===0) || (repO==="T1" && J===1) || (repO==="T2" && J===2) || (repO==="E" && J===2),
op=OSub[rep,r,J,JD,J13,TG,ND];
If[!MemberQ[opList,op] && op =!= 0,AppendTo[opList,op]]]],{J,Abs[JG-ND],JG+ND}],{TG,{"S","P","V","A"}}],{ND,0,MaxND}],
(* Moving frame *)
Do[Do[
JG=If[TG==="S" || TG==="P",0,If[TG==="V" || TG==="A",1,-100]];
Do[Do[
JD=0;
J13=0;
op=OSub[group,rep,r,J, \[Lambda],JD,J13,TG,ND];
If[!MemberQ[opList,op] && op =!= 0,AppendTo[opList,op]],{\[Lambda],0,J}],{J,Abs[JG-ND],JG+ND}],{TG,{"S","P","V","A"}}],{ND,0,MaxND}]];
(* Abs[JG-ND] could change to 0 if there are two derivatives *)
Return[opList];
];


OneHadronOperator[rep_,r_,MaxND_]:=OneHadronOperator[{0,0,0},rep,r,MaxND];


End[];
