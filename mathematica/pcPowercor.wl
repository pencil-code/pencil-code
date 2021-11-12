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
       Need not start from t=0. Uses FFT.
Output:
  {{0,1},{t2,ac2},{t3,ac3},...}"

fitTime::usage="fitTime[ts,lmodel,nfit] fits a time series using its first nfit points and 
a model specified by lmodel.
Input:
  ts:  List. Time series.
  lmodel:  String. Currently availabel models: \"EXP\", \"EXPCOS\", \"HALF\", \"AUTO\".
  nfit:  Integer.
Output:
  The decaying time as a real number."

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


autoCor[ts_]:=Module[{t,f,fft,ac},
  {t,f}=Transpose[ts];
  f=Transpose[Transpose[f]-Mean[f]];
  fft=If[Depth[f]==2,Fourier[f],Transpose[Fourier/@Transpose[f]]];
  ac=Re@InverseFourier[Total/@(Conjugate[fft]*fft)];
  {t-t[[1]],ac/ac[[1]]}//Transpose//Take[#,Length[t]/2//Floor]&
]

autoCorEnsemble[ts_,npiece_Integer:8]:=autoCor/@Partition[ts,Floor[Length[ts]/npiece]]


(* ::Section:: *)
(*Fitting*)


fitTime[ts_List,"EXP",-1]:=Module[{a,x},
  a[1]/.FindFit[ts//valueNorm,1*Exp[-x/a[1]],{a[1],a[2]},x]
]
fitTime[ts_List,"EXPCOS",-1]:=Module[{a,x},
  a[1]/.FindFit[ts//valueNorm,Cos[a[2]*x]*Exp[-x/a[1]],{a[1],a[2]},x]
]
fitTime[ts_List,"EXPCOSOMEGA",-1]:=Module[{a,x},
  a[2]/.FindFit[ts//valueNorm,Cos[a[2]*x]*Exp[-x/a[1]],{a[1],a[2]},x]
]
fitTime[ts_List,"HALF",nFit_Integer]:=Module[{t,v,pos,t1,t2,v1,v2},
  {t,v}=Transpose[Take[ts,nFit]];
  pos=Nearest[v->"Index",0.5,2];
  If[Abs[Subtract@@pos]!=1,t[[pos//Min]]//Return];
  {t1,t2}=Extract[t,List/@pos];{v1,v2}=Extract[v,List/@pos];
  (t1-t2+2t2*v1-2t1*v2)/(2v1-2v2)
]
fitTime[ts_List,"FFTSIN",nFit_Integer]:=Module[{t,v,fft},
  {t,v}=Transpose[ts];
  fft=Table[{i,Sin[i*t] . v},{i,0,2\[Pi]/Mean[Differences[t]],2\[Pi]/Last[t]}];
  2\[Pi]/6/MaximalBy[fft,Last][[1,1]]
]
fitTime[ts_List,"AUTO",nFit_Integer]:=Cases[
  fitTime[ts,#,nFit]&/@{"EXP","EXPCOS","HALF","FFTSIN"},
  x_?Positive
]//Min

(*using different models for k< or >kf*)
fitTime[ts_List,"EXPKF",nFit_Integer,knorm_:1,kf_:0]:=
  If[knorm<=kf+0.5,fitTime[ts,"EXPCOS",nFit],fitTime[ts,"EXP",nFit]
]
fitTime[ts_List,lFit_,nFit_Integer,knorm_,kf_]:=fitTime[ts,lFit,nFit]

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
    selectGood[tss_]:=If[nFit==-1,True,
      If[Union[Take[tss,nFit]//Differences//Transpose//Last//Sign]=={-1},True,False]];
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
  autoCor,fitTime,
  showResetTime,
  corrTime,corrTimeAC
]


EndPackage[]
