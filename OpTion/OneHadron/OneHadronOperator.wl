(* ::Package:: *)

OneHadronOperator::usage = "[ptot,rep,r,J,\[Lambda],JD,J13,TG,ND] generates one-hadron operators.";
OneHadronOperatorAll::usage = "OneHadronOperatorAll[ptot,rep,r,MaxND] generates all possible one-hadron operators.";


Begin["`OneHadronOperator`"];


OneHadronOperator[ptot_,rep_,r_,J_,\[Lambda]_,JD_,J13_,TG_,ND_]:=Module[{ope,group},
group=MomToGroup[ptot];
If[group==="Oh",
ope=OSub[rep,r,J,JD,J13,TG,ND],
ope=OSub[group,rep,r,J,\[Lambda],JD,J13,TG,ND]
];
Return[ope];
];


(* Moving frame: Search all possible nonzero operators *)
(* Only for MaxND <= 1 for the current verison *)
OneHadronOperatorAll[ptot_,rep_,r_,MaxND_]:=Module[{group,JG,JD,J13,repparity,repO,P,PGD,opList,op},
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
