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


BeginPackage["pc`","pcUtils`",
  "pcReadBasic`","pcRead1D`","pcRead2D`","pcReadFFT`",
  "pcReadVAR`","pcGetParam`","pcPlot`","pcPlot3D`",
  "pcTFM`","pcParticleStalker`","pcDerivative`",
  "pcPowercor`","pcPowerSpec`"
]


(* ::Chapter:: *)
(*Usage messages*)


pcInitialize::usage="pcInitialize[] is equivalent to (pcReload[];pcPlotStyle[]).";
pcParallelize::usage="pcParallelize[n] sets up the package for multiple kernels.
Input:
  n: Integer. Optional. Launches n subkernels. If not provided then launches all
     configured subkernels.
";

pcFunctions::usage="pcFunctions[] returns a list of all available functions in this package."


Begin["`Private`"]


(* ::Chapter:: *)
(*Functions*)


(* ::Section:: *)
(*One-click initialization*)


pcDirectory=$InputFileName//DirectoryName;
pcInitialize[]:=Module[{},
  pcReload[];pcPlotStyle[];
]
pcParallelize[n_|PatternSequence[]]:=With[{dir=pcDirectory},
  CloseKernels[];
  LaunchKernels[n]//Quiet;
  ParallelEvaluate[
    AppendTo[$Path,dir];
    Needs["pc`"];
    pcInitialize[]
  ];
]


(* ::Section:: *)
(*Extract all available functions*)


pcFunctions[]:=Cases[Names["pc*`*"],x_/;!StringContainsQ[x,"Private"]]


(* ::Chapter:: *)
(*End*)


End[]


Protect[
  pcInitialize,pcParallelize,
  pcFunctions
]


EndPackage[]
