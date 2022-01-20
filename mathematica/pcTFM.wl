(* ::Package:: *)

(* :Name: pcTFM` *)
(* :Author: Hongzhe Zhou, Stockholm, 2021*)
(* :Version: 0.1 *)
(* :Summary:
    Test-field method related functions.
*)


BeginPackage["pcTFM`","pcReadBasic`","pcRead1D`","pcUtils`"]


(* ::Chapter:: *)
(*Usage messages*)


TFMfindRange::usage="TFMfindRange[t,var,tReset,{shift,length},tsat,plotOption] can be
used for finding the appropriate range in the time series for extracting TFM results.
Input:
  t: time read from time_series.dat
  var: variable read from time_series.dat
  tReset: daainit
  {shift,length}: two Integers that specify the needed range
  tsat: only data after tsat is extracted
  plotOption: List, optional
Output:
  A ListPlot showing time series in black and selected data in red";

TFMResult::usage="TFMResult[sim,{varNames},{shift,length},tsat] gives a time series of the 
calculated transport coefficient.
Input:
  sim: String. Directory of the run
  varNames: Strings, one or more variables to be read from time_series.dat
  {shift,length}: two Integers that specify the needed range
  tsat: only data after tsat is extracted
Output:
  A time series of calculated transport coefficients";

TFMResultMean::usage="TFMResultMean[sim,{varNames},{shift,length},tsat] gives the mean value  
of the calculated transport coefficient.
Input:
  sim: String. Directory of the run
  varNames: Strings, one or more variables to be read from time_series.dat
  {shift,length}: two Integers that specify the needed range
  tsat: only data after tsat is extracted
Output:
  A number if data points is <37, and Around[mean,sd] otherwise";


Begin["`Private`"]


(* ::Chapter:: *)
(*Functions*)


(* ::Section:: *)
(*TFM*)


TFMextractTS[t_,var_,tReset_,{shift_,length_},tsat_]:=Module[
  {posReset,posPicked,posStart,posEnd,tsPicked},
  TFMextractTS::badparam="Chosen range too large. Returning {}.";
  posReset=1+First/@(Nearest[t->"Index",#]&/@Range[0.,t[[-1]],tReset]);
  posPicked=posReset/.x_Integer:>{x+shift,x+shift+length};
  If[Or@@Flatten[MapThread[Thread[#1>=#2]&,{Most@posPicked,Rest@posReset}]],
    Message[TFMextractTS::badparam];Return[{}]];
  posStart=First@FirstPosition[posPicked,x_List/;t[[x[[1]]]]>tsat];
  posEnd=First@FirstPosition[posPicked,x_List/;t[[x[[1]]]]>t[[-1]]-2tReset];
  posPicked=Take[posPicked,{posStart,posEnd}];
  Take[Transpose[{t,var}],#]&/@posPicked
]

TFMfindRange[t_,var_,tReset_,{shift_,length_},tsat_,plotOption___]:=ListLogPlot[
  Abs[{Transpose[{t,var}],Flatten[TFMextractTS[t,var,tReset,{shift,length},tsat],1]}]/.{x___,{}}:>{x},
  plotOption,ImageSize->Large,Joined->{True,False},PlotMarkers->{None,{Automatic,5}},
  PlotStyle->{Black,Red},GridLines->{{{tsat,{Red,Dashed,Thick}}},None}
]
 
TFMResult[sim_,{varNames__},{shift_,length_},tsat_]:=
  Module[{t,vars,tReset,coeffTS},
    t=readTS[sim,"t"];
    vars=readTS[sim,varNames];
    tReset=readParamNml[sim,"run.in","DAAINIT"];
    coeffTS[var_]:=pcAround/@TFMextractTS[t,var,tReset,{shift,length},tsat];
    coeffTS/@If[Length[Dimensions@vars]==1,{vars},vars]
  ]

TFMResultMean[sim_,{varNames__},{shift_,length_},tsat_]:=
  With[{vars=Map[Last[Transpose@#]&,TFMResult[sim,{varNames},{shift,length},tsat]]},
    Table[
      If[Length[data]<37,Mean[data],Around[Mean[#],StandardDeviation[#]]&@Map[Mean,Partition[data,6]]],
      {data,vars}
    ]
  ]


(* ::Chapter:: *)
(*End*)


End[]


Protect[
  TFMfindRange,TFMResult,TFMResultMean
]


EndPackage[]
