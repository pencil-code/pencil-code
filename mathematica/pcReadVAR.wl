(* ::Package:: *)

(* :Name: pcReadVAR` *)
(* :Author: Hongzhe Zhou, Stockholm, 2021*)
(* :Version: 0.1 *)
(* :Summary:
    Read VAR files.
*)


BeginPackage["pcReadVAR`","pcReadBasic`","pcDerivative`"]


(* ::Chapter:: *)
(*Usage messages*)


readVARN::usage="readVARN[sim,iVAR,addons,\"ltrim\"->True] reads the iVARth VARN file
from all processors.
Input:
  sim: String. Directory of the run
  iVAR: index of the VARN file, starting from 0
  addons: List. Optional. Specifies variables need to be computed; e.g., {\"oo\",\"bb\"}
  \"ltrim\": Options. Default value =True. If =True then trim ghost zones.
Output:
  An Association. Use readVARN[sim,iVAR]//Keys to extract its keys"

tSnap::usage="tSnap[sim,iproc:0] returns time of VARN files in the ith processor."

readVARNRaw;
readVARNProc;
dimProc;


Begin["`Private`"]


(* ::Chapter:: *)
(*Functions*)


(* ::Section:: *)
(*Parameter files for a single processor*)


dimProc[sim_,iproc_]:=Module[{mx,my,mz,c,precision,gh1,gh2,gh3,nx,ny,nz},
  {{mx,my,mz,c[1],c[2],c[3]},{precision},{gh1,gh2,gh3},c[4]}=Import@StringJoin[sim,"/data/proc",ToString@iproc,"/dim.dat"];
  precision=Switch[precision,"S","Real32","D","Real64"];
  {nx,ny,nz}=MapThread[(#1-2#2)&,{{mx,my,mz},{gh1,gh2,gh3}}];

  Thread[
    {"mx","my","mz","precision","gh1","gh2","gh3","nx","ny","nz"}->
    {mx,my,mz,precision,gh1,gh2,gh3,nx,ny,nz}
  ]//Association
]


gridProc[sim_,iproc_]:=Module[
  {file,mx,my,mz,precision,t,gridx,gridy,gridz,
  dx,dy,dz,Lx,Ly,Lz,dx1,dy1,dz1,dxt,dyt,dzt},
  gridProc::unfinished="Something left unread from `1`.";
  file=StringJoin[sim,"/data/proc",ToString@iproc,"/grid.dat"];
  {mx,my,mz,precision}=dimProc[sim,iproc]/@{"mx","my","mz","precision"};
  
  Close[file]//Quiet;
  BinaryRead[file,"Real32"];
  t=BinaryRead[file,precision];
  gridx=BinaryRead[file,ConstantArray[precision,mx]];
  gridy=BinaryRead[file,ConstantArray[precision,my]];
  gridz=BinaryRead[file,ConstantArray[precision,mz]];
  {dx,dy,dz}=BinaryRead[file,ConstantArray[precision,3]];
  BinaryRead[file,"Real32"];
  
  BinaryRead[file,"Real32"];
  {dx,dy,dz}=BinaryRead[file,ConstantArray[precision,3]];
  BinaryRead[file,"Real32"];
  
  BinaryRead[file,"Real32"];
  {Lx,Ly,Lz}=BinaryRead[file,ConstantArray[precision,3]];
  BinaryRead[file,"Real32"];
  
  BinaryRead[file,"Real32"];
  dx1=BinaryRead[file,ConstantArray[precision,mx]];
  dy1=BinaryRead[file,ConstantArray[precision,my]];
  dz1=BinaryRead[file,ConstantArray[precision,mz]];
  BinaryRead[file,"Real32"];
  
  BinaryRead[file,"Real32"];
  dxt=BinaryRead[file,ConstantArray[precision,mx]];
  dyt=BinaryRead[file,ConstantArray[precision,my]];
  dzt=BinaryRead[file,ConstantArray[precision,mz]];
  BinaryRead[file,"Real32"];
  
  If[BinaryRead[file,"Real32"]!=EndOfFile,Message[gridProc::unfinished,file]];
  Close[file];

  Flatten[Table[{gridx[[i]],gridy[[j]],gridz[[k]]},{k,mz},{j,my},{i,mx}],2]
]


tSnap[sim_,iproc_:0]:=
  With[{file=StringJoin[sim,"/data/proc",ToString@iproc,"/varN.list"]//Import},
    If[Head[file]===String,
      StringDrop[file,4]//ToExpression//List,
      file//Transpose//Last
    ]
  ]


(* ::Section:: *)
(*Read a single VARN (on a processor) file*)


readVARNProc[sim_,iproc_,iVAR_,nVar_,{pre_,mx_,my_,mz_},lshear_]:=
  Module[{file,var,typeLength,tmp,
    farray,t,gridx,gridy,gridz,dx,dy,dz,deltay,
    x,y,z},
        
    file=StringJoin[sim,"/data/proc",ToString@iproc,"/VAR",ToString@iVAR];
    Close[file]//Quiet;    
    var=BinaryReadList[OpenRead[file,BinaryFormat->True],Flatten@{
      "Integer32",ConstantArray[pre,mx*my*mz*nVar],"Integer32",
      (*assuming lshear=T; if not, deltay will be reset to 0 later*)
      "Integer32",ConstantArray[pre,1+mx+my+mz+3+1],"Integer32"
      },1,ByteOrdering->-1]//Flatten;
    (*warning: Something left unread, but we proceed*)
    Close[file];
    
    typeLength={1,mx*my*mz*nVar,1,1,1,mx,my,mz,1,1,1,1,1};
    {tmp,farray,tmp,tmp,t,gridx,gridy,gridz,dx,dy,dz,deltay,tmp}=Map[Take[var,#]&,
      MapThread[{#1-#2+1,#1}&,{typeLength//Accumulate,typeLength}]
    ];
    farray=Partition[farray,mx*my*mz];
    If[!lshear,deltay=0];
    
    {z,y,x}=Transpose@Flatten[Outer[List,gridz,gridy,gridx],2];
    Join[
      Thread[Keys[varName[sim]]->farray[[varName[sim]//Values]]],
      Thread[{"t","dx","dy","dz","deltay"}->Flatten@{t,dx,dy,dz,deltay}],
      Thread[{"gridx","gridy","gridz"}->{gridx,gridy,gridz}],
      Thread[{"x","y","z"}->{x,y,z}]
    ]//Association
  ]


(* ::Section:: *)
(*Combine a VARN file in all processors*)


readVARNRaw[sim_,iVAR_]:=With[{
  nproc=nProc[sim],
  nVar=Length@varName[sim],
  allvars=Join[varName[sim]//Keys,{"t","gridx","gridy","gridz","dx","dy","dz","deltay","x","y","z"}],
  lshear=readParamNml[sim,"start.in","lshear"]
  },
  Module[{data},
    data=Monitor[Table[
        iproc->readVARNProc[sim,iproc,iVAR,nVar,dimProc[sim,iproc]/@{"precision","mx","my","mz"},lshear],
        {iproc,0,nproc-1}
      ]//Association,
      StringJoin["Reading chunk ",ToString@iproc,"/",ToString@nproc]
    ];

  Table[var->Flatten[data[#][var]&/@Range[0,nproc-1]],{var,allvars}]//Association
]]

Options[readVARN]={"ltrim"->True}
readVARN[sim_,iVAR_,addOn_List:{},OptionsPattern[]]:=With[{
  vars=varName[sim]//Keys,
  raw=readVARNRaw[sim,iVAR],
  ltrim=OptionValue["ltrim"]
  },
  Module[{grid,values,merge,tmp,mx,my,mz,dx,dy,dz,gh1,gh2,gh3,trim,oo,bb,jj},
    readVARN::derivWOPeri="Warning: ghost zone values are incorrect when lperi!={T,T,T}";
    PrintTemporary["Combining files..."];
    
    grid=Transpose[raw/@{"x","y","z"}];
    values=raw/@vars;
    oo=bb=jj={{},{},{}};
    
    (*tmp=Transpose@SortBy[Transpose[{grid,Sequence@@values}]//DeleteDuplicates,First];*)
    (*tmp={grid,Sequence@@values}//Transpose//Union//Transpose;*)
    merge[{l__List}]/;Length[{l}]>=2:=Flatten[MaximalBy[#,Abs,1]&/@Transpose[{l}],1];
    merge[{l_List}]/;Length[{l}]==1:=l;
    tmp=Transpose[merge/@GatherBy[{grid,Sequence@@values}//Transpose//Union,First]];
    grid=tmp//First//Transpose;
    values=tmp//Rest;
    {mx,my,mz,gh1,gh2,gh3}=readDim[sim]/@{"mx","my","mz","gh1","gh2","gh3"};

    If[addOn!={},{dx,dy,dz}=Flatten[Union/@raw/@{"dx","dy","dz"}]];
    If[MemberQ[addOn,"oo"],
      oo=curl[values[[varName[sim]/@{"uu1","uu2","uu3"}]],mx,my,mz,gh1,gh2,gh3,dx,dy,dz]
    ];
    If[MemberQ[addOn,"bb"],
      bb=curl[values[[varName[sim]/@{"aa1","aa2","aa3"}]],mx,my,mz,gh1,gh2,gh3,dx,dy,dz]
    ];
    If[MemberQ[addOn,"jj"],
      jj=curl2[values[[varName[sim]/@{"aa1","aa2","aa3"}]],mx,my,mz,gh1,gh2,gh3,dx,dy,dz]
    ];
    
    If[ltrim,
      trim[f_]:=If[f=={},{},
        ArrayReshape[f,{mx,my,mz}][[gh1+1;;mx-gh1,gh2+1;;my-gh2,gh3+1;;mz-gh3]]//Flatten
      ];
      {grid,values,oo,bb,jj}=Map[trim,{grid,values,oo,bb,jj},{2}],
      (*else*)
      Message[readVARN::derivWOPeri]
    ];
      
    Join[
      Thread[{"t","dx","dy","dz","deltay"}->
        First/@(raw/@{"t","dx","dy","dz","deltay"})],
      Thread[{"lx","ly","lz"}->
        readDim[sim]/@If[ltrim,{"nx","ny","nz"},{"mx","my","mz"}]],
      Thread[{"x","y","z"}->grid],
      Thread[vars->values],
      (*long name to avoid shadowing the written auxillaries*)
      Thread[{"ooo1","ooo2","ooo3"}->oo],
      Thread[{"bbb1","bbb2","bbb3"}->bb],
      Thread[{"jjj1","jjj2","jjj3"}->jj]
    ]//Association
]]


(* ::Chapter:: *)
(*End*)


End[]


Protect[
  readVARN,tSnap
]


EndPackage[]
