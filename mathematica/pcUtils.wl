(* ::Package:: *)

(* :Name: pcUtils` *)
(* :Author: Hongzhe Zhou, Stockholm, 2022*)
(* :Version: 0.1 *)
(* :Summary:
    Here includes utils that are used by other parts of the package.
*)


BeginPackage["pcUtils`"]


(* ::Chapter:: *)
(*Usage messages*)


(*messages*)
pcAround::usage="pcAround[l] gives Around[l] if l={l1,l2,...},
or {Around[x],Around[y],...} if l={{x1,y1,...},{x2,y2,...},...}."

pcDivide::usage="pcDivide[l,n:3] partitions l into n pieces, with the
length of each no larger than Length[l]/n. The rest elements in l are
thrown away."

pcDifferences::usage="pcDifferences[l] returns
{(l2-l1)/2, (l2-l1)/2+(l4-l3)/2, ..., (l[n]-l[n-1])/2}."


Begin["`Private`"]


(* ::Chapter:: *)
(*Functions*)


pcAround[l_]:=Around[l]/;Depth[l]==2
pcAround[l_]:=Around/@Transpose[l]/;Depth[l]==3

pcDivide[l_,n_Integer:3]:=Partition[l,UpTo@Floor[Length[l]/n]][[1;;n]]

pcDifferences[l_]:=Module[{dl},
  dl=l//Differences;
  {dl[[1]]/2,Mean/@(Partition[dl,2,1]),dl[[-1]]/2}//Flatten
]


(* ::Chapter:: *)
(*End*)


End[]


Protect[
  pcAround,pcDivide,pcDifferences
]


EndPackage[]
