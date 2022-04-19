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

pcFit::usage="pcFit[data,sp,fact:1] fits the data with some model specified by
sp, prints out the result, and returns a the fitted curve.
Inputs:
  data: A List of the form {{x1,y1},{x2,y2},...,{xn,yn}.
  sp: Either a String that matches the implemented models, or simply an expression
      using \"x\" as the variable and \"a\"[i] as parameters.
      Example: \"a\"[1]+Log[\"a\"[2],\"x\"]^(\"a\"[3]).
      For the moment allows for up to 10 paramters.
  fact: The fitted curve is scaled by a factor of fact.
Outputs:
  A List of the form {{a1,b1},{a2,b2},...,{a32,b32}}, where a1=x1 and a32=xn."


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

pcFit[data_,sp_,fact_:1]:=Module[{model,a,x,sol,minmax},
  model=Switch[sp,
    "PowerLaw",a[1]*x^a[2],
    "Linear",a[1]+a[2]*x,
    "Exp",a[1]*Exp[a[2]*x],
    "Exp+C",a[1]+a[2]*Exp[a[3]*x],
    _,sp/.{"x"->x,"a"->a}
  ];
  sol=FindFit[data,model,a/@Range[10],x];
  Print["Fit result: ",model/.x->"x"/.sol];
  minmax=data[[;;,1]]//MinMax;
  Table[{x,fact*model/.sol},{x,Subdivide[Sequence@@minmax,32]}]
]


(* ::Chapter:: *)
(*End*)


End[]


Protect[
  pcAround,pcDivide,pcDifferences,pcFit
]


EndPackage[]
