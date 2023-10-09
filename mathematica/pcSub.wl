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
Also, use MollweideProjection[\"RegionFunction\",r,{thmin,thmax}] for the region
function."


Begin["`Private`"]


(* ::Chapter:: *)
(*Functions*)


MollweideProjection[theta_,phi_]:=Module[{lat,lon,aux,xx},
  lat=Pi/2-theta; lon=phi-Pi;
  aux=FindRoot[2xx+Sin[2xx]==Pi*Sin[lat],{xx,0.}][[1,2]];
  {2Sqrt[2]/Pi*lon*Cos[aux],Sqrt[2]*Sin[aux]}
]
MollweideProjection["RegionFunction",r0_,{thmin_,thmax_}]:=Function[{x,y,z},
  With[{th=ArcSin[y/r0/Sqrt[2]]},
    Between[ArcSin[(2th+Sin[2th])/Pi],{Pi/2-thmax,Pi/2-thmin}] &&
    Between[Pi*x/2/r0/Sqrt[2]/Cos[th],{-Pi,Pi}]
]]


(* ::Chapter:: *)
(*End*)


End[]


Protect[
  MollweideProjection
]


EndPackage[]
