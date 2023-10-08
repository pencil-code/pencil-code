(* ::Package:: *)

(* :Name: pcSub` *)
(* :Author: Hongzhe Zhou, Shanghai, 2023*)
(* :Version: 0.1 *)
(* :Summary:
    Mostly mathematical operations.
*)
(* :To-do:
    Move derivatives here.
*)


BeginPackage["pcSub`"]


(* ::Chapter:: *)
(*Usage messages*)


(*messages*)
MollweideProjection::usage="MollweideProjection[\[Theta],\[Phi]] gives the
Mollweide projection coordinate {x,y}, where \[Theta]\[Element][0,\[Pi]] and \[Phi]\[Element][0,2\[Pi]].
Also, use MollweideProjection[\"RegionFunction\",r] for the region
function."


Begin["`Private`"]


(* ::Chapter:: *)
(*Functions*)


MollweideProjection[theta_,phi_]:=Module[{aux,xx},
  aux=FindRoot[2xx+Sin[2xx]==Pi*Sin[Pi/2-theta],{xx,0.}][[1,2]];
  {2Sqrt[2]/Pi*(phi-Pi)*Cos[aux],Sqrt[2]*Sin[aux]}
]
MollweideProjection["RegionFunction",r0_]:=Function[{x,y,z},
  With[{th=ArcSin[y/r0/Sqrt[2]]},
    -Pi/2<=ArcSin[(2th+Sin[2th])/Pi]<=Pi/2 &&
    0<=Pi+Pi*x/2/r0/Sqrt[2]/Cos[th]<=2Pi
]]


(* ::Chapter:: *)
(*End*)


End[]


Protect[
  MollweideProjection
]


EndPackage[]
