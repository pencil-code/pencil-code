(* ::Package:: *)

(* :Name: pcPowercor` *)
(* :Author: Hongzhe Zhou, Stockholm, 2021*)
(* :Version: 0.1 *)
(* :Summary:
    Read data and compute correlation time from subroutine powercor.
*)


BeginPackage["pcPowercor`","pcRead1D`"]


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

corrTimePowercor::usage="corrTime[sim,tsat,file,lFit,nFit] returns the correlation times.
Input:
  sim: String. Run directory.
  tsat: NumericQ. Starting point of the saturated phase.
  file: String. powercor_*.dat or correlation_*.dat.
  lFit: String. Fitting model.
  nFit: Integer or All. Number of points to be fitted.
Output:
  {correlation time of the mean of all AC's, a List of correlation times for each AC}"


Begin["`Private`"]


(* ::Chapter:: *)
(*Functions*)


(* ::Section:: *)
(*Utils*)


(*shift an autocorrelation horizontally to start from t=0*)
timeShift[l_List]:=Transpose[Transpose[l]-{l[[1,1]],0}]
(*shift an autocorrelation vertically to start from value 1*)
valueNorm[l_List]:=Transpose[Transpose[l]/{1,l[[1,2]]}]
(*do both timeShift and valueNorm*)
normalizeAC[l_List]:=l//timeShift//valueNorm


(* ::Section:: *)
(*Auto-correlation*)


autoCor[ts_]:=Module[{t,f,fft,ac},
  {t,f}=Transpose[ts];
  f=Transpose[Transpose[f]-Mean[f]];
  fft=If[Depth[f]==2,Fourier[f],Transpose[Fourier/@Transpose[f]]];
  ac=Re@InverseFourier[Total/@(Conjugate[fft]*fft)];
  {t-t[[1]],ac}//Transpose//Take[#,Length[t]/2//Floor]&
]
autoCor[t_,var_]:=autoCor[{t,var}//Transpose]


(* ::Section:: *)
(*Fitting*)


fitTime[ts_List,model_String,nFit_]:=Module[{a,x},
  a[1]/.FindFit[Take[ts,nFit]//normalizeAC,
    Switch[model,
      "EXP",Exp[-x/a[1]],
      "EXP+C",1-a[2]+a[2]*Exp[-x/a[1]],
      "EXPCOS",Cos[a[2]*x]*Exp[-x/a[1]],
      "EXPCOSOMEGA",Cos[a[1]*x]*Exp[-x/a[2]]
    ],{a[1],a[2]},x]
]
fitTime[ts_List,"EXPPOS",nFit_]:=Module[{t,v,pos,x,a},
  {t,v}=Transpose[ts//normalizeAC];
  pos=FirstPosition[v,_?Negative,All];
  If[pos=!=All,pos=pos[[1]]-1;If[pos<=2,Return[-1]]];
  a[1]/.FindFit[Take[ts//normalizeAC,pos],Exp[-x/a[1]],a[1],x]
]
fitTime[ts_List,"HALF",nFit_]:=Module[{t,v,pos,t1,t2,v1,v2},
  {t,v}=Transpose[Take[ts,nFit]//normalizeAC];
  pos=Nearest[v->"Index",0.5,2];
  If[Abs[Subtract@@pos]!=1,t[[pos//Min]]//Return];
  {t1,t2}=Extract[t,List/@pos];{v1,v2}=Extract[v,List/@pos];
  (t1-t2+2t2*v1-2t1*v2)/(2v1-2v2)
]
fitTime[ts_List,"FFTSIN",nFit_]:=Module[{t,v,fft},
  {t,v}=Transpose[Take[ts,nFit]//normalizeAC];
  fft=Table[{i,Sin[i*t] . v},{i,0,2\[Pi]/Mean[Differences[t]],2\[Pi]/Last[t]}];
  2\[Pi]/6/MaximalBy[fft,Last][[1,1]]
]
fitTime[ts_List,"INT",nFit_]:=Module[{t,v},
  {t,v}=Transpose[Take[ts,nFit]//normalizeAC];
  Total[v]*Mean[Differences[t]]
]
fitTime[ts_List,"MIN",nFit_Integer]:=Cases[
  fitTime[ts,#,nFit]&/@{"EXP","EXPCOS","HALF","FFTSIN","INT"},
  x_?Positive
]//Min

(*using different models for k< or >kf*)
fitTime[ts_List,"EXPKF",nFit_,knorm_:1,kf_:0]:=
  If[knorm<=kf+0.5,fitTime[ts,"EXPCOS",nFit],fitTime[ts,"EXP",nFit]
]
fitTime[ts_List,lFit_,nFit_,knorm_,kf_]:=fitTime[ts,lFit,nFit]


(* ::Section:: *)
(*Find correlation time*)


corrTimePowercor[sim_,tsat_,file_,lFit_,nFit_]:=Module[{tvart,t,spec,pos,d},
  (*find reset time of sp*)
  tvart=Import[sim<>"/data/tvart.dat"]//Flatten//Cases[#,tt_/;tt>tsat]&//Most;
  (*read spectra*)
  {t,spec}=read1D[sim,file];
  (*find positions of tvart*)
  pos=Nearest[t->"Index",#]&/@tvart//Flatten;
  (*cut spec into intervals*)
  d=pos//Differences//Min;
  pos={pos,pos+d-1}//Transpose;
  spec=Transpose@Mean[Take[spec,#]&/@pos];
  (*form time series*)
  spec=Transpose[{Range[0,d-1]*(t//Differences//Mean),#}]&/@spec;
  (*fitting*)
  {fitTime[spec//Mean,lFit,nFit],fitTime[#,lFit,nFit]&/@spec}
]


(* ::Chapter:: *)
(*End*)


End[]


Protect[
  autoCor,fitTime,
  corrTimePowercor
]


EndPackage[]
