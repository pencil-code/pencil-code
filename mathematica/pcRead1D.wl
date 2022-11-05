(* ::Package:: *)

(* :Name: pcRead1D` *)
(* :Author: Hongzhe Zhou, Stockholm, 2021*)
(* :Version: 0.1 *)
(* :Summary:
    Read time series and other 1D data files.
*)


BeginPackage["pcRead1D`","pcReadBasic`"]


(* ::Chapter:: *)
(*Usage messages*)


readTS::usage="readTS[sim,Options] reads time_series.dat and returns an Association object.
It always returns rectangular data by filling the missing elements by Missing[] (e.g. in
cases where print.in has been modified before continuing a run).
Input:
  sim: String. Directory of the run
Options:
  \"lRest\": By default True, which removes the first record in the file (usually when t=0)
Output:
  An Association object whose keys include all the variables that have appeared
Examples:
  ts=readTS[sim], and then ts[\"t\"]={t1,t2,...}, ts[\"urms\"]={u1,u2,...}, etc.
  readTS[sim,var] is the same as readTS[sim][var].
  readTS[sim,var1,var2,...] is the same as readTS[sim]/@{var1,var2,...}";

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
  file: String. Name of the data file, including .dat.
        Can also be multiple Strings; then the file name will be str1<>str2<>... .
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

readSpecFluxPQ::usage="readSpecFluxPQ[sim] returns the wavenumber coordinates for the
spectral flux files, based on the values of specflux_dp and specflux_dq in run.in."

readSpecFlux::usage="readSpecFlux[sim,file] reads spectral flux files.
Input:
  sim: String. Directory of the run
  file: String. Name of the data file, including .dat
Output:
  {{t1, t2, ..., tn}, d}, where d is a n\[Cross]npq\[Cross]3 List."


Begin["`Private`"]


(* ::Chapter:: *)
(*Functions*)


(* ::Section:: *)
(*Time series*)


readTSSingle[ts_List]:=Module[{head,value},
  head=Rest@StringSplit[ts[[1]],"-"..];
  value=ts//Rest//StringReplace[#,{"e"->"*^","E"->"*^"}]&//StringSplit//ToExpression;
  AssociationThread[head->Transpose[value]]
]

readTS[sim_String,OptionsPattern[{"lRest"->True}]]:=Module[{str,ts,fullHead,lookUp},
  (* read time_series.dat into a 1D List of Strings *)
  str=OpenRead@FileNameJoin[{sim,"data","time_series.dat"}];
  ts=ReadList[str,Record];
  Close[str];
  (* for each chunk startring with "#", reform it into an Association *)
  ts=readTSSingle/@Split[ts,!StringStartsQ[#2,"#"]&];
  
  (* merge chunks and fill with Missing[]*)
  fullHead=Union@Flatten[Keys/@ts];
  lookUp[key_]:=Lookup[#,key,ConstantArray[Missing[],#[[1]]//Length]]&/@ts;
  
  If[OptionValue["lRest"],
    AssociationMap[Rest@*Flatten@*lookUp,fullHead],
    AssociationMap[Flatten@*lookUp,fullHead]
  ]
]
readTS[sim_String,var_String,opt:OptionsPattern[]]:=readTS[sim,opt][var]
readTS[sim_String,vars__String,opt:OptionsPattern[]]:=readTS[sim,opt]/@{vars}

growthRate[sim_,var_,t1_,t2_]:=Module[{t,f,pos,a,tt},
{t,f}=readTS[sim,"t",var];
pos=Nearest[t->"Index",#]&/@{t1,t2}//Flatten;
{t,f}=Take[#,pos]&/@{t,f};
a[1]/.FindFit[Transpose[{t,Log[f]}],a[1]*tt+a[2],{a[1],a[2]},tt]
]


(* ::Section:: *)
(*Spectra-like files*)


read1D[sim_String,file__String]:=Module[{str,ts},
  str=OpenRead[FileNameJoin[{sim,"data",StringJoin[file]}]];
  ts=ReadList[str,Number,RecordLists->True];
  Close[str];
  
  ts=Flatten/@Split[ts,Length[#2]>1&];
  
  {ts[[;;,1]],ts[[;;,2;;]]}
]
read1D[sim_String,file__String,l_Integer]:=Module[{str,ts},
  str=OpenRead[FileNameJoin[{sim,"data",StringJoin[file]}]];
  ts=ReadList[str,Number,RecordLists->False]//Partition[#,l+1]&;
  Close[str];
  
  {ts[[;;,1]],ts[[;;,2;;]]}
]


read1DSigned[sim_,file_,l_Integer:0,testF_:Positive]:=Module[{t,spec},
  {t,spec}=If[l<=0,read1D[sim,file],read1D[sim,file,l]];
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


(* ::Section:: *)
(*Spectral fluxes*)


readSpecFluxPQ[sim_]:=Module[{pmin,pmax,dp,dq,nk,grid},
  pmin=readParamNml[sim,"run.in","SPECFLUX_PMIN"];
  pmax=readParamNml[sim,"run.in","SPECFLUX_PMAX"];
  dp=readParamNml[sim,"run.in","SPECFLUX_DP"];
  dq=readParamNml[sim,"run.in","SPECFLUX_DQ"];
  nk=readDim[sim]["nx"]/2;
  
  grid[d_,min_:0,max_:nk-1]:=If[d>0,Range[min,max,d],(-d)^Range[0,Floor@Log[-d,max]]];
  
  Round@{grid[dp,pmin,pmax],grid[dq]}
]
readSpecFlux[sim_,file_]:=Module[{p,q,t,tra,tmp},
  {p,q}=readSpecFluxPQ[sim];
  {t,tra}=read1D[sim,file,Length[p]*Length[q]];
  Do[
    tmp[i][p[[j]],q[[k]]]=Partition[tra[[i]],Length[p]][[k,j]],
    {k,q//Length},{j,p//Length},{i,t//Length}];
  
  {t,p,q,Table[tmp[i],{i,Length[t]}]}
]


(* ::Chapter:: *)
(*End*)


End[]


Protect[
  readTS,growthRate,
  read1D,read1DSigned,read1D2Scale,
  readAves,
  readSpecFluxPQ,readSpecFlux
]


EndPackage[]
