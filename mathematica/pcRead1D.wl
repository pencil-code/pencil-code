(* ::Package:: *)

(* :Name: pcRead1D` *)
(* :Author: Hongzhe Zhou, Stockholm, 2021*)
(* :Version: 0.1 *)
(* :Summary:
    Read time series and other 1D data files.
*)


BeginPackage["pcRead1D`"]


(* ::Chapter:: *)
(*Usage messages*)


readTS::usage="readTS[sim,var] reads the column in time_series.dat that corresponds to 
variable var. Support multiple vars: readTS[dir,var1,var2,...].
Input:
  sim: String. Directory of the simulation folder
  var: String. The variable to read; e.g., \"t\", or \"urms\".
       var=\"HEAD\" gives all possible entries
Output:
  { {t1,t2,...}, {var1,...}, {var2,...} }";

growthRate::usage="growthRage[sim,var,t1,t2] computes the exponential growth rate
of var between time t1 and t2.
Input:
  sim: String. Directory of the simulation folder
  var: String. The variable to read
  t1,t2: NumericQ
Output:
  A real number"

read1D::usage="read1D[sim,file,l] reads 1D data.
Input:
  sim: String. Directory of the simulation folder
  file: String. Name of the data file, including .dat
  l:  Integer. Optional. If exists then specifies the length of the 1D data.
Output:
  {{t1,t2,...},{{data at t1},{data at t2},...}}"

read1DSigned::usage="read1DSigned[sim,file,l:0,testF:Positive] reads 1D data and
separates them into groups based on the test function testF. Those with testF
returning False are replaced by zeros in the first group, and those with testF
returning True are replaced by zeros in the second group.
Input:
  sim: String. Directory of the simulation folder
  file: String. Name of the data file, including .dat
  l:  Integer. Optional. If exists then specifies the length of the 1D data.
  testF: The test function. testF==Positive by default. Other examples: Real, #^2<0.5&
Example output with testF==Positive:
  { {t1,t2,...},
    {{1,0,2,3,0,...},...},
    {{0,-4,0,0,-5,...},...}
  }"

read1D2Scale::usage="read1D2Scale[sim,file,kf] reads 1D spectrum and splits it into
large- and small-scale parts according to the dividing wave number kf.
Input:
  sim: String. Directory of the simulation folder
  file: String. Name of the data file, including .dat
  kf: Integer. Modes with wave numbers >= kf will be included in the small-scale part.
Output:
  {t, large-scale spectra, small-scale spectra}.";

readAves::usage="readAves[sim,sp,varNames] returns planar averaged values.
Input:
  sim: String. Directory of the run
  sp: String. \"xy\", \"yz\", or \"xz\", specifying which plane is averaged over
  varNames: Strings, one or more variables to be read. Need to be registered in sim/xyaver.in
Output:
  A List with its first element {t1,t2,...}, and the rest time series of Nx(yz) dimentional data."


Begin["`Private`"]


(* ::Chapter:: *)
(*Functions*)


(* ::Section:: *)
(*Time series*)


readTS[sim_,vars__]:=Module[{file,ts,head,pos},
  readTS::badform="`1`: Inconsistent number of columns in ts.";
  readTS::novar="`1`: Variable `2` not found in ts.";
  file=sim<>"/data/time_series.dat";
  ts=DeleteCases[Import[file],x_/;Length[x]==1];
  If[Not[Equal@@(Length/@ts)],
    Message[readTS::badform,sim];
    Return["BadTsFormat"]
  ];
  head=Rest@Flatten[StringSplit[First@Import[file],"-"..]];
  If[vars=="HEAD",Return[head]];
  If[Length[ts[[1]]]!=Length[head],
    Message[readTS::badform,sim];
    Return[$Failed]
  ];
  Table[
    If[(pos=Position[head,var])=={},
      Message[readTS::novar,sim,var];
      ConstantArray[$Failed,ts//Length],
      (*else*)
      ts[[;;,pos[[1,1]]]]
    ],
    {var,{vars}}
  ]/.{{x__}}:>{x}
]

growthRate[sim_,var_,t1_,t2_]:=Module[{t,f,pos,a,tt},
{t,f}=readTS[sim,"t",var];
pos=Nearest[t->"Index",#]&/@{t1,t2}//Flatten;
{t,f}=Take[#,pos]&/@{t,f};
a[1]/.FindFit[Transpose[{t,Log[f]}],a[1]*tt+a[2],{a[1],a[2]},tt]
]


(* ::Section:: *)
(*Spectra-like files*)


read1D[sim_,file_]:=Module[{data,pos},
  data=Import[sim<>"/data/"<>file];
  pos=Flatten@Position[data,x_/;Length[x]==1];
  data=Flatten/@Partition[data,pos[[2]]-pos[[1]]];
  {First/@data,Rest/@data}
]
read1D[sim_,file_,l_]:=Module[{data},
  data=Import[sim<>"/data/"<>file]//Flatten;
  data=Partition[data,l+1];
  {First/@data,Rest/@data}
]

read1DSigned[sim_,file_,l_Integer:0,testF_:Positive]:=Module[{t,spec},
  {t,spec}=If[l<=0,read1D[sim,file],read1D[sim,file,l]];
  spec=Rest/@spec;
  {t,spec/.x_?(!testF[#]&)->0,spec/.x_?testF->0}
]

Options[read1D2Scale]={"lSum"->False};
read1D2Scale[sim_,file_,kf_Integer:-1,OptionsPattern[]]:=Module[{t,spec,n,specL,specS,kff},
  {t,spec}=read1D[sim,file];
  kff=Round@If[kf==-1,
    If[Length[#[[1]]]==1,#[[2,1]],#[[1,2]]]&@Import[FileNameJoin[{sim,"k.dat"}]],
    kf
  ];
  n=spec//First//Length;
  specL=(UnitStep[kff-#]&/@Range[n])*#&/@spec;
  specS=(UnitStep[#-kff-1]&/@Range[n])*#&/@spec;
  If[OptionValue["lSum"],
    {t,Total/@specL,Total/@specS},
    {t,specL,specS}
  ]
]


(* ::Section:: *)
(*Planar averaged files*)


readAves[sim_,plane_,varNames__]:=
  With[{inFile=sim<>"/"<>plane<>"aver.in",datFile=sim<>"/data/"<>plane<>"averages.dat"},
    Module[{varList,index,t,vars,error=False},
      readAves::nofile="No .in or .dat file found. Returning {0,0}.";
      readAves::novar="Variable `1` not found. Returning {0,0}.";
      If[!And[FileExistsQ@inFile,FileExistsQ@datFile],error=True;Message[readAves::nofile]];
      varList=Flatten[{Import@inFile}];
      Scan[If[!MemberQ[varList,#],error=True;Message[readAves::novar,#]]&,Flatten@{varNames}];
      If[error,Return[{0,0}]];
      
      index=Flatten@Map[Position[varList,#]&,Flatten@{varNames}];
      {t,vars}=read1D[sim,plane<>"averages.dat"];
      vars=Transpose@Map[Partition[#,Length[#]/Length[varList]]&,vars];
      Join[{t},vars[[index]]]
    ]
  ]


(* ::Chapter:: *)
(*End*)


End[]


Protect[
  readTS,growthRate,
  read1D,read1DSigned,read1D2Scale,
  readAves
]


EndPackage[]
