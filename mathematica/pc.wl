(* ::Package:: *)

(* :Name: pc` *)
(* :Author: Hongzhe Zhou, Stockholm, 2021*)
(* :Version: 0.1 *)
(* :Summary:
    Needs this package automatically Needs the following:
      pcReadBasic, pcRead1D, pcRead2D,
      pcReadVAR, pcGetParam, pcPlot,
      pcTFM, pcParticleStalker,pcDerivative,
      pcPowercor
*)


BeginPackage["pc`",
  "pcReadBasic`","pcRead1D`","pcRead2D`",
  "pcReadVAR`","pcGetParam`","pcPlot`",
  "pcTFM`","pcParticleStalker`","pcDerivative`",
  "pcPowercor`"
]


(* ::Chapter:: *)
(*Usage messages*)


pcInitialize::usage="pcInitialize[] is equivalent to (pcReload[];pcPlotStyle[]).";

pcFunctions::usage="pcFunctions[] returns a list of all available functions in this package."


Begin["`Private`"]


(* ::Chapter:: *)
(*Functions*)


(* ::Section:: *)
(*One-click initialization*)


pcInitialize[]:=(pcReload[];pcPlotStyle[];)


(* ::Section:: *)
(*Extract all available functions*)


pcFunctions[]:=Cases[Names["pc*`*"],x_/;!StringContainsQ[x,"Private"]]


(* ::Chapter:: *)
(*End*)


End[]


Protect[
  pcInitialize,
  pcFunctions
]


EndPackage[]
