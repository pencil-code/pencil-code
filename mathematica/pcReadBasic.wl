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
  
readDim::usage="readDim[sim,n] reads dimension-related parameters.
Input:
  sim: String. Directory of the simulation folder
  n: Integer. Optional. If absent then read the global file; otherwise read procn/dim.dat.
Output:
  An Association object with available arguments
    {mx,my,mz,precision,gh1,gh2,gh3,nx,ny,nz}.";

readGrid::usage="readGrid[sim] reads from /data/allprocs/grid.dat only for the xyz grids.
Input:
  sim: String. Directory of the simulation folder
Output:
  {gridx,gridy,gridz}, each being a 1D List.";

varName::usage="varName[sim] returns contents in data/varName.dat.
Input:
  sim: String. Directory of the simulation folder
Output:
  ";

nProc::usage="nProc[sim] returns a list of number of processors in xyz directions."
procBounds::usage="procBounds[sim] returns proc bounds."

nSnap::usage="nSnap[sim,iproc:0] returns to number of VARN files in the ith processor."


Begin["`Private`"]


(* ::Chapter:: *)
(*Functions*)


(* ::Section:: *)
(*Read basic information*)


readDim[sim_,n_Integer:-1]:=Module[{file,data,precision,nxyz,names},
  file=StringJoin[sim,"/data",If[n==-1,"","/proc"<>ToString@n],"/dim.dat"];
  data=Import[file];
  precision=Switch[data[[2,1]],"S","Real32","D","Real64"];
  nxyz=MapThread[(#1-2#2)&,{data[[1,1;;3]],data[[3]]}];
  
  names=Join[{"mx","my","mz","mvar","maux","mglobal","precision","gh1","gh2","gh3"},
    If[n==-1,{"nprocx","nprocy","nprocz","Iprocz_slowest"},{"ipx","ipy","ipz"}],
    {"nx","ny","nz"}];
  
  AssociationThread[names->Flatten@List[data[[1]],precision,data[[3;;-1]],nxyz]]
]

Options[readGrid]={"ltrim"->True};
readGrid[sim_,iproc_Integer:-1,OptionsPattern[]]:=Module[{file,p,mx,my,mz,tmp},
  file=If[iproc==-1,
    StringJoin[sim,"/data/allprocs/grid.dat"],
    StringJoin[sim,"/data/proc"<>ToString[iproc]<>"/grid.dat"]
  ];
  {p,mx,my,mz}=readDim[sim,iproc]/@{"precision","mx","my","mz"};
  
  file//Close//Quiet;
  BinaryRead[file,"Integer32"];
  BinaryRead[file,p]; (* time *)
  tmp=BinaryRead[file,ConstantArray[p,#]]&/@{mx,my,mz};
  Close[file];
  
  If[OptionValue["ltrim"],tmp[[;;,4;;-4]],tmp]
]

nProc[sim_]:=Import[sim<>"/data/dim.dat"]//Last//Most

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
    "T",  True,
    "F",  False,
    str_/;StringContainsQ[str,LetterCharacter],value,
    str_/;StringStartsQ[str,DigitCharacter],ToExpression[value],
    str_/;StringStartsQ[str,"+"],ToExpression[value],
    str_/;StringStartsQ[str,"-"],ToExpression[value],
    _,    value
  ]
]


(* ::Chapter:: *)
(*End*)


End[]


Protect[
  readParamNml,readDim,readGrid,
  varName,
  nProc,procBounds,
  nSnap
]


EndPackage[]
