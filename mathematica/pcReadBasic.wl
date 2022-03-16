(* ::Package:: *)

(* :Name: pcReadBasic` *)
(* :Author: Hongzhe Zhou, Stockholm, 2021*)
(* :Version: 0.1 *)
(* :Summary:
    Functions that read files.
*)


BeginPackage["pcReadBasic`"]


(* ::Chapter:: *)
(*Usage messages*)


readParamNml::usage="readParamNml[sim,file,var] reads the value of the parameter p in the file f
var from param namelist files.
Input:
  sim: String. Directory of the simulation folder
  file: String. Either \"start.in\" or \"run.in\"
  var: String. Need to match an entry in file; e.g., \"NU\"
Output:
  A value, or String flags like T, F, NaN, or others";
  
readDim::usage="readDim[sim] reads dimension-related parameters.
Input:
  sim: String. Directory of the simulation folder
Output:
  An Association object with available arguments
    {mx,my,mz,precision,gh1,gh2,gh3,nx,ny,nz}.";

varName::usage="varName[sim] returns contents in data/varName.dat.
Input:
  sim: String. Directory of the simulation folder
Output:
  ";

nProc::usage="nProc[sim] returns the total number of processors."
procBounds::usage="procBounds[sim] returns proc bounds."

nSnap::usage="nSnap[sim,iproc:0] returns to number of VARN files in the ith processor."


Begin["`Private`"]


(* ::Chapter:: *)
(*Functions*)


(* ::Section:: *)
(*Read basic information*)


readDim[sim_]:=Module[{mx,my,mz,c,precision,gh1,gh2,gh3,nx,ny,nz},
  {{mx,my,mz,c[1],c[2],c[3]},{precision},{gh1,gh2,gh3},c[4]}=Import@StringJoin[sim,"/data/dim.dat"];
  precision=Switch[precision,"S","Real32","D","Real64"];
  {nx,ny,nz}=MapThread[(#1-2#2)&,{{mx,my,mz},{gh1,gh2,gh3}}];
  
  Thread[
    {"mx","my","mz","precision","gh1","gh2","gh3","nx","ny","nz"}->
    {mx,my,mz,precision,gh1,gh2,gh3,nx,ny,nz}
  ]//Association
]

nProc[sim_]:=Length@FileNames[sim<>"/data/proc*"]-1

procBounds[sim_]:=Module[{pre,preByte,file,nx,ny,nz,bx,by,bz},
  pre=readDim[sim]["precision"];
  preByte=If[pre=="Real32",4,8];
  file=sim<>"/data/proc_bounds.dat";
  Close[file]//Quiet;
  nx=BinaryRead[file,"Integer32"]/preByte;
  bx=BinaryRead[file,ConstantArray[pre,nx]];
  BinaryRead[file,"Integer32"];
  ny=BinaryRead[file,"Integer32"]/preByte;
  by=BinaryRead[file,ConstantArray[pre,ny]];
  BinaryRead[file,"Integer32"];
  nz=BinaryRead[file,"Integer32"]/preByte;
  bz=BinaryRead[file,ConstantArray[pre,nz]];
  BinaryRead[file,"Integer32"];
  Close[file];
  
  {bx,by,bz}
]

nSnap[sim_,iproc_:0]:=Import[sim<>"/data/tsnap.dat"][[1,2]]

varName[sim_]:=Rule@@@Reverse/@Import[sim<>"/data/varname.dat"]//Association


(* ::Section:: *)
(*Read simulation parameters*)


readParamNml[sim_,file_,param_]:=Module[{breakTwoInOne,listToString,import,value},
  readParamNml::noentry="Entry `1` not found for `2`.";
  
  (*import nml file; one line normally goes to one List; if one line becomes a String, break it*)
  breakTwoInOne[nml_]:=Replace[Import[sim<>"/data/"<>nml],x_String:>StringSplit[x,","],{1}];
  (*for each List in x, make it a long String*)
  listToString[x_]:=StringReplace[StringRiffle[#,","]," "..->""]&/@x;
  import[nml_]:=listToString[breakTwoInOne[nml]];

  value=Cases[
    Switch[file,
      "start.in",import["param.nml"],
      "run.in",import["param2.nml"]
    ],_?(StringMatchQ[#,param<>"=*"]&)
  ];
  
  If[value=={},Message[readParamNml::noentry,param,sim];Return[$Failed]];
  value=StringReplace[First[value],x__~~",":>x];
  If[StringContainsQ[value,","],Return[
      StringSplit[StringDrop[value,StringLength[param]+1],","]
    ]
  ];
  
  value=StringDrop[
    StringReplace[value,{"E+"->"*^+","E-"->"*^-"}],StringLength[param]+1
  ];
  Switch[value,
    str_/;StringStartsQ[str,DigitCharacter],ToExpression[value],
    str_/;StringStartsQ[str,"+"],ToExpression[value],
    str_/;StringStartsQ[str,"-"],ToExpression[value],
    "T",  True,
    "F",  False,
    _,    value
  ]
]


(* ::Chapter:: *)
(*End*)


End[]


Protect[
  readParamNml,readDim,
  varName,
  nProc,procBounds,
  nSnap
]


EndPackage[]
