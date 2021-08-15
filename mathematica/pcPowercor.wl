(* ::Package:: *)

(* :Name: pcPowercor` *)
(* :Author: Hongzhe Zhou, Stockholm, 2021*)
(* :Version: 0.1 *)
(* :Summary:
    Read data and compute correlation time from subroutine powercor.
*)


BeginPackage["pcPowercor`"]


(* ::Chapter:: *)
(*Usage messages*)


autoCor::usage="autoCor[ts] computes auto-correlation of ts.
Input:
  ts:  List. A time series data; for example readTS[sim,\"t\",\"urms\"]//Transpose.
       Need not start from t=0.
Output:
  {{0,1},{t2,ac2},{t3,ac3},...}"

autoCorEnsemble::usage="autoCorEnsemble[ts,n:8] cuts ts into 8 pieces and then
computes the auto-correlation of each of them.
Input:
  ts:  List. A time series data; for example readTS[sim,\"t\",\"urms\"]//Transpose.
       Need not start from t=0.
  n:  Size of the ensemble. By default n=8.
Output:
  A List of length n, each of the form {{0,1},{t2,ac2},{t3,ac3},...}"

showResetTime::usage="showResetTime[{t,f},dtcor,shift:0,plotStyle:{}] returns a ListPlot that
shows both the original curve and the identified resetting points.
Input:
  t,f:  List
  dtcor:  Resetting time
  shift:  Integer. Optional. Shift red dots.
  plotStyle:  List. Optional
Output:
  A ListPlot where the original curve is black and reseting points are red."

corrTime::usage="corrTime[{t,f},dtcor,model,nFit,\"lSingle\"->False] returns the correlation time.
If t and f are long enough then produces a mean with its standard deviation.
For now neighboring 5 curves are binned and the size of the ensemble is at least 3.
Input:
  t,f:  List. Data produced with dtcor
  dtcor:  Resetting time
  model:  String. Specifies model: \"EXP\", \"EXPCOS\", \"LIN\"
  nFit:  Integer. Number of points to take to fit
  lSingle:  Optional. If True then forces to produce a single value
Output:
  Around[mean,sd] or mean depending on the length of t and f"

corrTimeAC::usage="corrTime[{t,f},model,nFit,npiece:8] returns the correlation
time computed from the auto-correlation function formed by {t,f}.
Input:
  t,f:  List. Time series data.
  model:  String. Specifies model: \"EXP\", \"EXPCOS\", \"LIN\"
  nFit:  Integer. Number of points to take to fit
  npiece:  Size of the ensemble. If =1 then produces a single value
Output:
  Around[mean,sd] or mean depending on npiece"


Begin["`Private`"]


(* ::Chapter:: *)
(*Functions*)


(* ::Section:: *)
(*Utils*)


mean[x_List]:=With[{y=Cases[x,e_?NumericQ]},
  If[y=={},Return["Bad curve"]];
  If[Length[y]==1,y[[1]],Around[Mean[y],StandardDeviation[y]]]
]

timeShift[l_List]:=With[{t0=l[[1,1]]},l/.{x_,y_}:>{x-t0,y}]
valueNorm[l_List]:=With[{l0=l[[1,2]]},l/.{x_,y_}:>{x,y/l0}]
signChange[l_List]:=With[{s=Sign[l[[1,2]]]},l/.{x_,y_}:>{x,y/s}]

(*findResetTime[t_,dtcor_]:=Rest[Most[
  Nearest[t->"Index",0.05+#,2][[2]]&/@Range[0,t[[-1]],dtcor]
]]*)

findResetTime[t_,dtcor_]:=Module[{pos},
  pos=Position[t,tt_/;tt>dtcor && Mod[Round[tt],Round[dtcor]]==1]//Flatten;
  Extract[pos,Position[Differences[pos],d_/;d>5]]
]

showResetTime[{t_,f_},dtcor_,shift_Integer:0,plotStyle_List:{}]:=With[
  {pos=findResetTime[t,dtcor],ts=Transpose[{t,f}]},
  ListPlot[{ts,ts[[pos+shift]]},Join[plotStyle,{
      ImageSize->800,Joined->{True,False},
      PlotStyle->{Black,Red},PlotRange->All
      }]
  ]
]

splitCurve[t_,f_,dtcor_]:=With[{pos=findResetTime[t,dtcor]},Module[{pos1,t1,f1,lmin},
  pos1=Partition[pos,2,1]/.{x_?NumericQ,y_}:>{x,y-1};
  t1=Take[t,#]&/@pos1;
  f1=Take[f,#]&/@pos1;
  lmin=Max[16,Min[Length/@f1]];
  timeShift/@MapThread[Take[Transpose[{#1,#2}],lmin]&,{t1,f1}]
]]


(* ::Section:: *)
(*Auto-correlation*)


autocor[f_,\[Tau]_,dt_]:=With[{times=If[Length[Dimensions[f]]==1,Times,Dot]},
  MapThread[times[#1,#2]*#3&,{Drop[f,-\[Tau]],Drop[RotateLeft[f,\[Tau]],-\[Tau]],Drop[dt,-\[Tau]]}]//Total
]
autoCor[ts_]:=With[{ts0=timeShift[ts]},Module[{t,f,dt,norm,\[Tau]max},
  {t,f}=Transpose[ts0];
  f=Most[f];
  dt=Differences[t];
  \[Tau]max=Round[Length[dt]/3];
  norm=autocor[f,0,dt];
  {t[[1;;\[Tau]max+1]],autocor[f,#,dt]/norm&/@Range[0,\[Tau]max]}//Transpose
]]

autoCorEnsemble[ts_,npiece_Integer:8]:=autoCor/@Partition[ts,Floor[Length[ts]/npiece]]


(* ::Section:: *)
(*Fitting*)


fitTime[ts_List,model_String,nFit_Integer]:=Module[{a,x},
  a[1]/.FindFit[Take[ts,nFit]//valueNorm,
    Switch[model,
      "EXP",a[2]*Exp[-x/a[1]],
      "EXPCOS",a[2]*Cos[a[3]*x]*Exp[-x/a[1]]
    ],{a[1],a[2],a[3]},x]
]
fitTime[ts_String,model_String,nFit_Integer]:=ts


(* ::Section:: *)
(*Find correlation time*)


Options[corrTime]={"lSingle"->False,"lSelectGood"->True}
corrTime[{t_List,f_List},dtcor_,model_String,nFit_Integer,OptionsPattern[]]:=With[{
  nBin=5, (*bin neighboring 5 curves*)
  nMinCurve=3 (*min size of ensemble*)
  },Module[{selectGood,ts=splitCurve[t,f,dtcor]}, 
    selectGood[tss_]:=If[Union[Take[tss,nFit]//Differences//Transpose//Last//Sign]=={-1},True,False];
    If[Length[ts[[-1]]]<nFit,ts=Most[ts]];
    ts=signChange/@ts;
    If[OptionValue["lSelectGood"],ts=Cases[ts,x_?selectGood]];
    If[ts=={},Return["Bad curve."]];
    If[Length[ts]<nBin*nMinCurve || OptionValue["lSingle"],
      fitTime[Mean[ts],model,nFit],
      mean[fitTime[Mean[#],model,nFit]&/@Partition[ts,nBin]]
    ]
]]

corrTimeAC[{t_List,f_List},model_String,nFit_Integer,npiece_Integer:8]:=
  With[{ts=Transpose[{t,f}]},
  If[npiece==1,
    fitTime[autoCor[ts],model,nFit],
    mean[fitTime[#,model,nFit]&/@autoCorEnsemble[ts,npiece]]
  ]
]


(* ::Chapter:: *)
(*End*)


End[]


Protect[
  autoCor,autoCorEnsemble,
  showResetTime,
  corrTime,corrTimeAC
]


EndPackage[]
