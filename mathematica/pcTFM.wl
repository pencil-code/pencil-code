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


TFMfindRange::usage="TFMfindRange[sim,tsat,var,{shift,length},plotOption] can be
used for finding the appropriate range in the time series for extracting TFM results.
Input:
  sim: Run directory of the simulation.
  tsat: only data after tsat is extracted.
  var: variable read from time_series.dat
  {shift,length}: two Integers that specify the needed range
  plotOption: Optional. Additional plot style options.
Output:
  A ListPlot showing time series in black and selected data in red";

TFMResult::usage="TFMResult[sim,tsat,vars,{shift,length}] gives a time series of the 
calculated transport coefficient.
Input:
  sim: String. Directory of the run
  tsat: only data after tsat is extracted
  vars: Can be one String, or a List of Strings.
  {shift,length}: two Integers that specify the needed range
Output:
  A time series of calculated transport coefficients";

TFMResultMean::usage="TFMResultMean[sim,tsat,{vars},{shift,length}] gives the mean value  
of the calculated transport coefficient.
Input:
  sim: String. Directory of the run
  tsat: only data after tsat is extracted
  vars: Can be one String, or a List of Strings.
  {shift,length}: two Integers that specify the needed range
Output:
  A number if data points is <37, and Around[mean,sd] otherwise";

TFMResultCmplx::usage="TFMResultCmplx[sim,tsat,var,{shift,length}] returns the real and
imaginary parts of a coefficient, given that OM_TESTFIELD!=0.
Input:
  sim: String. Directory of the run
  tsat: only data after tsat is extracted
  var: String, one variable to be read from time_series.dat
  {shift,length}: two Integers that specify the needed range
Output:
  {Re[var],Im[var]}";


Begin["`Private`"]


(* ::Chapter:: *)
(*Functions*)


(* ::Section:: *)
(*TFM*)


getDINIT[sim_,var_]:=With[{sp=StringTake[var,3]},Which[
  MemberQ[{"alp","eta","b11","b12","b21","b22"},sp],"DAAINIT",
  MemberQ[{"gam","kap","c1r","c2r"},sp],"DCCINIT"
  ]//readParamNml[sim,"run.in",#]&
]

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

TFMfindRange[sim_,tsat_,var_,{shift_,length_},plotOption___]:=Module[{t,f,dinit},
  {t,f}=readTS[sim,"t",var];
  dinit=getDINIT[sim,var];
  ListPlot[
    ({Transpose[{t,f}],Flatten[TFMextractTS[t,f,dinit,{shift,length},tsat],1]})/.{x___,{}}:>{x},
    plotOption,ImageSize->Large,Joined->{True,False},PlotMarkers->{None,{Automatic,5}},
    PlotStyle->{Black,Red},GridLines->{{{tsat,{Red,Dashed,Thick}}},None},
    ScalingFunctions->"SignedLog"
  ]
]

TFMResult[sim_,tsat_,var_String,{shift_,length_}]:=Map[
  pcAround,
  TFMextractTS[readTS[sim,"t"],readTS[sim,var],getDINIT[sim,var],{shift,length},tsat]
]
TFMResult[sim_,tsat_,vars_List,{shift_,length_}]:=Map[
  pcAround/@TFMextractTS[readTS[sim,"t"],readTS[sim,#],getDINIT[sim,#],{shift,length},tsat]&,
  vars
]

TFMResultMean[sim_,tsat_,var_String,{shift_,length_}]:=Mean[Last/@TFMResult[sim,tsat,var,{shift,length}]]
TFMResultMean[sim_,tsat_,vars_List,{shift_,length_}]:=Map[
  Mean,
  Map[Last,TFMResult[sim,tsat,vars,{shift,length}],{2}]
]

TFMResultCmplx[sim_,tsat_,var_String,{shift_,length_}]:=Module[{t,x,dt,om},
  {t,x}=TFMResult[sim,tsat,var,{shift,length}]//Transpose;
  dt=t//Differences//Mean;
  om=readParamNml[sim,"run.in","OM_TESTFIELD"];
  Through[
    {Re,Im}[dt*(Exp[I*om*t] . x)/(t[[-1]]-t[[1]])]
  ]
]


(* ::Chapter:: *)
(*End*)


End[]


Protect[
  TFMfindRange,TFMResult,TFMResultMean,TFMResultCmplx
]


EndPackage[]
