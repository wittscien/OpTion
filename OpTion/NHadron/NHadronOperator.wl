(* ::Package:: *)

NHadronOperator::usage = "NHadronOperator[ptot,rep,r,momTuple,parTuple,msTuple] generates N-hadron operators by projection.";
NHadronOperatorAll::usage = "NHadronOperatorAll[ptot,rep,r,maxmom,parTuple] generates all possible N-hadron operators by projection.";


Begin["`NHadronOperator`"];


(* Combine all functions. For this N-Hadron version, it is just a wrap-around ProjectionOperatorN.  *)
NHadronOperator[ptot_,rep_,r_,momTuple_,parTuple_,msTuple_]:=Module[{group,op,Npar,sTuple,pTuple,parTupleOrder},
(* Sanity check *)
Npar=Length[momTuple];
Assert[N===Length[parTuple]];
Assert[N===Length[msTuple]];
Assert[ptot===Total[momTuple]];
(* Extract the particle information and reorganize *)
sTuple=ConstantArray[0,Npar];
pTuple=ConstantArray[0,Npar];
parTupleOrder=ConstantArray[0,Npar];
Do[If[Head[parTuple[[i]]]===List,{sTuple[[i]],pTuple[[i]]}=parTuple[[i]],{sTuple[[i]],pTuple[[i]]}=ParToSpinParity[parTuple[[i]]]];parTupleOrder[[i]]={sTuple[[i]],pTuple[[i]]},{i,Npar}];
(* Translates the total momentum into the (little) group *)
(* 2024.08.19: No more needed *)
(*group=If[IntegerQ[Total[sTuple]],MomToGroup[ptot],MomToGroup[ptot,1]];*)
(* Call the funtion *)
op=ProjectionOperatorN[ptot,rep,r,momTuple,parTupleOrder,msTuple];
Return[op];
];


(* Search for all operators *)
NHadronOperatorAll[ptot_,rep_,r_,maxmom_,parTuple_]:=Module[{group,opList,op,Npar,sTuple,pTuple,momList},
(* Extract the particle information and reorganize *)
Npar=Length[parTuple];
sTuple=ConstantArray[0,Npar];
pTuple=ConstantArray[0,Npar];
Do[If[Head[parTuple[[i]]]===List,{sTuple[[i]],pTuple[[i]]}=parTuple[[i]],{sTuple[[i]],pTuple[[i]]}=ParToSpinParity[parTuple[[i]]]],{i,Npar}];
group=If[IntegerQ[Total[sTuple]],MomToGroup[ptot],MomToGroup[ptot,1]];
momList=GenerateMomVectorsM[maxmom,Npar,ptot];
(* Search for all possible spins and momenta *)
opList={};
Do[
op=NHadronOperator[ptot,rep,r,momTuple,parTuple,msTuple];
If[!MemberQ[opList,op] && op =!= 0 && IsLinearlyIndependent[opList,op],AppendTo[opList,op]],{momTuple,momList},{msTuple,Tuples[Table[Range[-sTuple[[pari]],sTuple[[pari]]],{pari,Npar}]]}];
Return[opList];
];


NHadronOperatorAll[rep_,r_,maxmom_,parTuple_]:=NHadronOperatorAll[{0,0,0},rep,r,maxmom,parTuple];


End[];
