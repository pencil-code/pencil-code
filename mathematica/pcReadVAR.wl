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


readVARNProc::usage="readVARNProc[sim,iproc,iVAR,nVar,lVar:{},lshear,Options] reads a VAR file
from a single processor. There is now a new and faster implementation readVARProc.
Input:
  sim: String. Directory of the run.
  iproc: Index of the processor, starting from 0.
  iVAR: Integer. Index of the VARN file, starting from 0. If iVAR==-1 then will read var.dat,
        i.e., the last snapshot.
  nVAR: Usually Length[varName[sim]].
  lVar: Optional. By default {}. One can put in integers like varName[sim][\"lnrho\"].
        This will instruct to only read certain variables. For the moment does
        not support magic variables like \"oo1\" or \"bb1\".
  lshear: True for a shear run, and false for a non-shear run.
Options:
  \"lnoghost\": By default ->True. Will remove all the *inner* ghost zones. The ones at the
                outermost region of the simulation region will not be removed.
Output:
  An Association. Use Keys to extract its keys.
Example:
  readVARNProc[sim,1,2,7,{\"lnrho\"},False] will read the \"lnrho\" variable in proc1/VAR2 for
  a non-shear run. The inner ghost zones will be removed but not the outermost ones."

readVARProc::usage="readVARProc[sim,iproc,iVAR,globalInfo] is a new and faster implementation
of readVARNProc. It reads a VAR file from a single processor.
Input:
  sim: String. Directory of the run.
  iproc: Index of the processor, starting from 0.
  iVAR: Integer. Index of the VARN file, starting from 0. If iVAR==-1 then will read var.dat,
        i.e., the last snapshot.
  globalInfo: List. Optional. If provided it should be {nProc[sim],varName[sim]//Keys}.
Options:
  \"ltrimInner\": Trim the inner ghost zones or not. False by default.
  \"ltrimOuter\": Trim the outer ghost zones or not. False by default.
  \"imagic\": List of magic variables. Possible ones are \"oo\", \"bb\", and \"jj\".
              They will be labeled by upcase letters like \"OO1\",\"OO2\",\"OO3\" in the output.
Output:
  An Association object.
Example:
  readVARProc[sim,1,2,\"ltrimInner\"->True,\"ltrimOuter\"->True,\"imagic\"->{\"oo\",\"jj\"}]"

readVARN::usage="readVARN[sim,iVAR,addOn,\"ltrim\"->False,\"lOnly\"->False] reads the iVARth VARN file
from all processors. Works for io_dist. For io_hdf5, simply use Import[\".../VAR1.h5\"].
Input:
  sim: String. Directory of the run
  iVAR: Integer. Index of the VARN file, starting from 0. If iVAR==-1 then will read var.dat,
        i.e., the last snapshot.
  addOn: List. Optional. Specifies magical variables need to be computed; e.g., {\"oo\",\"bb\"}.
         Magical variables in the result are assigned with long names, e.g., \"ooo1\" and \"bbb1\"
         to avoid shadowing the written auxiliaries.
  \"ltrim\": Option. Default value ->False.
             If ->True then trim ghost zones.
             If ->\"More\" then trim again nghost number of cells at the boundaries; i.e., 2*nghost cells are
               removed at each boundary.
  \"lOnly\": Option. Default value ->Flase.
             If ->True, then will only read those variables in addOn.
Output:
  An Association. Use readVARN[sim,iVAR]//Keys to extract its keys
Examples:
  readVARN[sim,0] reads all variables in VAR0, and does not trim any ghost zones on the boundary.
  readVARN[sim,0,{\"bb\"},\"ltrim\"->True] not only read all variables, but also computes the three components
    of the magnetic field. Furthermore, ghost zones are trimmed.
  readVARN[sim,0,{\"bb\"},\"ltrim\"->\"More\",\"lOnly\"->True] will only return the three components of the
    magnetic field. Further more, 2*ghost zones are trimmed on the boundary. For the moment it is not possible
    to read only one component of the magnetic field. Similarly for oo and jj.
  readVARN[sim,0,{\"uu1\"},\"lOnly\"->True] will only read the uu1 component.
  readVAR[sim,0,{\"uu1\"},\"lOnly\"->False] is the same as readVAR[sim,0]."

readVAR::usage="readVAR[sim,iVAR] is a new and faster implementation of readVARN. It reads a VAR
file for all the processors. Works for io_dist. For io_hdf5, simply use Import[\".../VAR1.h5\"].
Input:
  sim: String. Directory of the run
  iVAR: Integer. Index of the VARN file, starting from 0. If iVAR==-1 then will read var.dat,
        i.e., the last snapshot.
Options:
  \"ltrim\": Trim the outer ghost zones or not. True by default.
  \"imagic\": List of magic variables. Possible ones are \"oo\", \"bb\", and \"jj\".
              They will be labeled by upcase letters like \"OO1\",\"OO2\",\"OO3\" in the output.
Output:
  An Association object.
Examples:
  readVAR[sim,0,\"ltrim\"->False,\"imagic\"->{\"oo\",\"jj\"}]"

tSnap::usage="tSnap[sim,iproc:0] returns time of VARN files in the ith processor."

pcIntpltVar::usage="pcIntpltVar[VAR,vars] gives interpolated functions of vars in the VAR file.
Input:
  VAR: The returned object of readVARN[...].
  vars: A List of Strings, e.g., {\"bbb1\",\"bbb2\",\"bbb3\"}.
Output:
  A List of InterpolatingFunctions."


Begin["`Private`"]


(* ::Chapter:: *)
(*Functions*)


(* ::Section:: *)
(*Parameter files for a single processor*)


gridProc[sim_,iproc_]:=Module[
  {file,mx,my,mz,precision,t,gridx,gridy,gridz,
  dx,dy,dz,Lx,Ly,Lz,dx1,dy1,dz1,dxt,dyt,dzt},
  gridProc::unfinished="Something left unread from `1`.";
  file=StringJoin[sim,"/data/proc",ToString@iproc,"/grid.dat"];
  {mx,my,mz,precision}=readDim[sim,iproc]/@{"mx","my","mz","precision"};
  
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
(*Read all variables in a single VARN (on a processor) file*)


Options[readVARNProc]={"lnoghost"->True};
readVARNProc[sim_String,iproc_Integer,iVAR_Integer,nVar_Integer,lVar_List:{},lshear_,OptionsPattern[]]:=
  Module[{file,stream,var,typeLength,tmp,varPos,
    pre,mx,my,mz,gh1,gh2,gh3,ipx,ipy,ipz,
    farray,t,gridx,gridy,gridz,dx,dy,dz,deltay,
    x,y,z,trim,posx,posy,posz,nproc},
    
    {pre,mx,my,mz,gh1,gh2,gh3,ipx,ipy,ipz}=readDim[sim,iproc]/@{"precision","mx","my","mz","gh1","gh2","gh3","ipx","ipy","ipz"};
    nproc=nProc[sim];
    file=If[iVAR==-1,
      StringJoin[sim,"/data/proc",ToString@iproc,"/var.dat"],
      StringJoin[sim,"/data/proc",ToString@iproc,"/VAR",ToString@iVAR]
    ];
    Close[file]//Quiet;
    stream=OpenRead[file,BinaryFormat->True];
    If[lVar=={},
      (*read all variables*)
      var=BinaryReadList[stream,Flatten@{
        "Integer32",ConstantArray[pre,mx*my*mz*nVar],"Integer32",
        "Integer32",ConstantArray[pre,1+mx+my+mz+3+1],"Integer32"
        },1,ByteOrdering->-1]//Flatten;
      Close[file];
      typeLength={1,mx*my*mz*nVar,1,1,1,mx,my,mz,1,1,1,1,1};
      {tmp,farray,tmp,tmp,t,gridx,gridy,gridz,dx,dy,dz,deltay,tmp}=
        Map[Take[var,#]&,MapThread[{#1-#2+1,#1}&,{typeLength//Accumulate,typeLength}]],
      (*else, only read required ones*)
      varPos=Map[4+ToExpression[StringTake[pre,-2]]/8*mx*my*mz*(#-1)&,lVar];
      var=Flatten@{
        (SetStreamPosition[stream,#];BinaryRead[stream,ConstantArray[pre,mx*my*mz]])&/@varPos,
        SetStreamPosition[stream,4+ToExpression[StringTake[pre,-2]]/8*mx*my*mz*nVar+4+4];
        BinaryRead[stream,ConstantArray[pre,1+mx+my+mz+3+1]]
      };
      Close[file];
      typeLength={mx*my*mz*Length@lVar,1,mx,my,mz,1,1,1,1};
      {farray,t,gridx,gridy,gridz,dx,dy,dz,deltay}=
        Map[Take[var,#]&,MapThread[{#1-#2+1,#1}&,{typeLength//Accumulate,typeLength}]]
    ];
    farray=Partition[farray,mx*my*mz];
    If[!lshear,deltay=0];
    
    (*remove inner ghost zones*)
    If[OptionValue["lnoghost"]!=True, trim[x_]:=x,
      posx=Which[
        ipx==0,1;;mx-gh1,
        ipx==nProc[sim][[1]]-1,gh1+1;;-1,
        True,gh1+1;;mx-gh1
      ];
      posy=Which[
        ipy==0,1;;my-gh2,
        ipy==nProc[sim][[2]]-1,gh2+1;;-1,
        True,gh2+1;;my-gh2
      ];
      posz=Which[
        ipz==0,1;;mz-gh3,
        ipz==nProc[sim][[3]]-1,gh3+1;;-1,
        True,gh3+1;;mz-gh3
      ];
      If[nproc[[1]]==1,posx=All];
      If[nproc[[2]]==1,posy=All];
      If[nproc[[3]]==1,posz=All];
      
      (*trim*)
      gridx=gridx[[posx]];
      gridy=gridy[[posy]];
      gridz=gridz[[posz]];
      trim[x_]:=ArrayReshape[x,{mz,my,mx}][[posz,posy,posx]]//Flatten
    ];
    
    (*the index is still f[[mz,my,mz]] here*)
    {z,y,x}=Transpose@Flatten[Outer[List,gridz,gridy,gridx],2];
    Join[
      Thread[Keys[varName[sim]][[If[lVar=={},All,lVar]]]->trim/@farray],
      Thread[{"t","dx","dy","dz","deltay"}->Flatten@{t,dx,dy,dz,deltay}],
      Thread[{"gridx","gridy","gridz"}->{gridx,gridy,gridz}],
      Thread[{"x","y","z"}->{x,y,z}]
    ]//Association
  ]


(* ::Section:: *)
(*Combine a VARN file in all processors*)


Options[readVARNRaw]={"lnoghost"->True};
readVARNRaw[sim_,iVAR_,lVar_List:{},OptionsPattern[]]:=Module[{
  nproc=Times@@nProc[sim],nVar=Length@varName[sim],lshear=readParamNml[sim,"start.in","lshear"],
  data,allvars},
  
  allvars=Join[
    Part[varName[sim]//Keys,If[lVar=={},All,lVar]],
    {"t","gridx","gridy","gridz","dx","dy","dz","deltay","x","y","z"}];
  data=Monitor[Table[
    iproc->readVARNProc[sim,iproc,iVAR,nVar,lVar,lshear,"lnoghost"->OptionValue["lnoghost"]],
    {iproc,0,nproc-1}]//Association,
    StringJoin["Reading chunk ",ToString@iproc,"/",ToString@nproc]
  ];

  Table[var->Flatten[data[#][var]&/@Range[0,nproc-1]],{var,allvars}]//Association
]

Options[readVARN]={"ltrim"->False,"lOnly"->False,"i2d"->Automatic};
readVARN[sim_,iVAR_,addOn_List:{},OptionsPattern[]]:=Module[{
  ltrim=OptionValue["ltrim"],lonly=OptionValue["lOnly"],
  lVar,vars=varName[sim]//Keys,varPos,raw,i2d,
  grid,values,tmp,mx,my,mz,dx,dy,dz,gh1,gh2,gh3,trim,trimx,trimy,trimz,
  oo,bb,jj},
  readVARN::newone="Warning: there is a now a faster and more robust implementation readVAR.";
  readVARN::ghvalues="Warning: ghost zone values are computed assuming periodic boundary conditions.";
  readVARN::trimmore="Error: Cannot do ltrim=More because of insufficient inner grid cells.";
  
  Message[readVARN::newone];
  
  (*If lOnly==True, we read only variables in addOn*)
  If[lonly,
    (*find non-magic variables*)
    lVar=DeleteCases[addOn,_?(StringMatchQ[{"oo*","bb*","jj*"}])];
    (*add necessary variables to lVar*)
    If[MemberQ[addOn,"oo"],lVar={lVar,"uu1","uu2","uu3"}];
    If[MemberQ[addOn,"bb"],lVar={lVar,"aa1","aa2","aa3"}];
    If[MemberQ[addOn,"jj"],lVar={lVar,"aa1","aa2","aa3"}];
    lVar=Sort@Flatten[Position[vars,#]&/@Union[lVar//Flatten]];
    raw=readVARNRaw[sim,iVAR,lVar,"lnoghost"->True],
    
    (*else, we read all*)
    lVar=All;raw=readVARNRaw[sim,iVAR,"lnoghost"->True]
  ];
  vars=vars[[lVar]];
  varPos=AssociationThread[vars->Range[Length@vars]];
  
  PrintTemporary["Combining files..."];
  grid=Transpose[raw/@{"x","y","z"}];
  values=raw/@vars;
  oo=bb=jj={{},{},{}};
  
  (*Sort by grid. Now indexing is f[[mx,my,mz]] for everything*)
  tmp={grid,Sequence@@values}//Transpose//SortBy[#,First]&//Transpose;
  grid=tmp//First//Transpose;
  values=tmp//Rest;
  {mx,my,mz,gh1,gh2,gh3}=readDim[sim]/@{"mx","my","mz","gh1","gh2","gh3"};
  
  If[addOn!={},{dx,dy,dz}=Flatten[Union/@raw/@{"dx","dy","dz"}]];
  If[MemberQ[addOn,"oo"],oo=curl[values[[varPos/@{"uu1","uu2","uu3"}]],mx,my,mz,gh1,gh2,gh3,dx,dy,dz]];
  If[MemberQ[addOn,"bb"],bb=curl[values[[varPos/@{"aa1","aa2","aa3"}]],mx,my,mz,gh1,gh2,gh3,dx,dy,dz]];
  If[MemberQ[addOn,"jj"],jj=curl2[values[[varPos/@{"aa1","aa2","aa3"}]],mx,my,mz,gh1,gh2,gh3,dx,dy,dz]];
  
  (*trim ghost zones*)
  i2d=If[OptionValue["i2d"]===Automatic,
    If[Times@@(readDim[sim]/@{"nprocx","nprocy","nprocz"})>1,-1,1],
    OptionValue["i2d"]
  ];
  trimx=If[mx==1+2*gh1,i2d,gh1+1;;mx-gh1];
  trimy=If[my==1+2*gh2,i2d,gh2+1;;my-gh2];
  trimz=If[mz==1+2*gh3,i2d,gh3+1;;mz-gh3];
  trim[{}]={};
  trim[f_]:=Switch[ltrim,
    False,f,
    True, ArrayReshape[f,{mx,my,mz}][[trimx,trimy,trimz]]//Flatten,
    "More",If[(mx-2gh1)*(my-2gh2)*(mz-2gh3)<=0,Message[readVARN::trimmore];Return[$Failed]];
           ArrayReshape[f,{mx,my,mz}][[2gh1+1;;mx-2gh1,2gh2+1;;my-2gh2,2gh3+1;;mz-2gh3]]//Flatten
  ];
  If[ltrim==False,Message[readVARN::ghvalues]];
  {grid,values,oo,bb,jj}=Map[trim,{grid,values,oo,bb,jj},{2}];
  
  tmp=Join[
    Thread[{"t","dx","dy","dz","deltay"}->First/@(raw/@{"t","dx","dy","dz","deltay"})],
    Thread[{"x","y","z"}->grid],
    Thread[vars->values]
  ]//Association;
  (*long name to avoid shadowing the written auxiliaries*)
  If[oo!={{},{},{}},tmp=Join[tmp,AssociationThread[{"ooo1","ooo2","ooo3"}->oo]]];
  If[bb!={{},{},{}},tmp=Join[tmp,AssociationThread[{"bbb1","bbb2","bbb3"}->bb]]];
  If[jj!={{},{},{}},tmp=Join[tmp,AssociationThread[{"jjj1","jjj2","jjj3"}->jj]]];
  tmp
]


(* ::Section:: *)
(*New implementations of readVARProc and readVAR*)


Options[readVARProc]={"ltrimInner"->False,"ltrimOuter"->False,"imagic"->{}};
readVARProc[sim_,iproc_,iVAR_,globalInfo_List:{},OptionsPattern[]]:=Module[{
  mx,my,mz,pre,ipx,ipy,ipz,gh1,gh2,gh3,
  nprocx,nprocy,nprocz,varList,nvar,dx,dy,dz,
  file,stream,read,var,t,
  magic,tmp,
  ltrimI,ltrimO
},
  
  With[{dim=readDim[sim,iproc]},
    {mx,my,mz,pre,gh1,gh2,gh3}=dim/@{"mx","my","mz","precision","gh1","gh2","gh3"};
    {ipx,ipy,ipz}=dim/@{"ipx","ipy","ipz"};
  ];
  
  {{nprocx,nprocy,nprocz},varList}=If[globalInfo=={},
    {nProc[sim],varName[sim]//Keys},
    globalInfo
  ];
  nvar=varList//Length;
  
  (* on a single proc, VAR contains: mz*mx*my data of nvar, t,x,y,z,dx,dy,dz,(deltay) *)
  file=If[iVAR==-1,
    StringJoin[sim,"/data/proc",ToString@iproc,"/var.dat"],
    StringJoin[sim,"/data/proc",ToString@iproc,"/VAR",ToString@iVAR]
  ];
  Close[file]//Quiet;
  stream=OpenRead[file,BinaryFormat->True];
  read[type_List]:=BinaryReadList[stream,type,1,ByteOrdering->-1][[1]];
  read[type_]:=BinaryRead[stream,type,ByteOrdering->-1];
  
  read["Integer32"];
  var=read[ConstantArray[pre,mx*my*mz*nvar]]//ArrayReshape[#,{nvar,mz,my,mx}]&;
  read[{"Integer32","Integer32"}];
  t=read[pre];
  read[ConstantArray[pre,mx+my+mz]];
  {dx,dy,dz}=read[{pre,pre,pre}];
  Close[file];
  
  var=Transpose[var,{1,4,3,2}];
  
  (* compute magic variables oo, bb, or jj *)
  magic=OptionValue["imagic"];
  If[MemberQ[magic,"oo"],
    varList=Flatten@{varList,"OO1","OO2","OO3"};
    tmp=curl[
      var[[varName[sim]/@{"uu1","uu2","uu3"}]]//ArrayReshape[#,{3,mx*my*mz}]&,
      mx,my,mz,gh1,gh2,gh3,dx,dy,dz
    ]//ArrayReshape[#,{3,mx,my,mz}]&;
    var=Join[var,tmp]
  ];
  If[MemberQ[magic,"bb"],
    varList=Flatten@{varList,"BB1","BB2","BB3"};
    tmp=curl[
      var[[varName[sim]/@{"aa1","aa2","aa3"}]]//ArrayReshape[#,{3,mx*my*mz}]&,
      mx,my,mz,gh1,gh2,gh3,dx,dy,dz
    ]//ArrayReshape[#,{3,mx,my,mz}]&;
    var=Join[var,tmp]
  ];
  If[MemberQ[magic,"jj"],
    varList=Flatten@{varList,"JJ1","JJ2","JJ3"};
    tmp=curl2[
      var[[varName[sim]/@{"aa1","aa2","aa3"}]]//ArrayReshape[#,{3,mx*my*mz}]&,
      mx,my,mz,gh1,gh2,gh3,dx,dy,dz
    ]//ArrayReshape[#,{3,mx,my,mz}]&;
    var=Join[var,tmp]
  ];
  
  (* take out ghost zones *)
  ltrimI=OptionValue["ltrimInner"];
  ltrimO=OptionValue["ltrimOuter"];
  If[(ltrimI && ipz<nprocz-1 ) || (ltrimO && ipz==nprocz-1), var=Map[Drop[#,-gh3]&,var,{1}]];
  If[(ltrimI && ipz>0 )        || (ltrimO && ipz==0),        var=Map[Drop[#,gh3]&,var,{1}]];
  If[(ltrimI && ipy<nprocy-1 ) || (ltrimO && ipy==nprocy-1), var=Map[Drop[#,-gh2]&,var,{2}]];
  If[(ltrimI && ipy>0 )        || (ltrimO && ipy==0),        var=Map[Drop[#,gh2]&,var,{2}]];
  If[(ltrimI && ipx<nprocx-1 ) || (ltrimO && ipx==nprocx-1), var=Map[Drop[#,-gh1]&,var,{3}]];
  If[(ltrimI && ipx>0 )        || (ltrimO && ipx==0),        var=Map[Drop[#,gh1]&,var,{3}]];
  
  Join[
    Association["ipxyz"->{ipx,ipy,ipz}],
    Association["t"->t],
    Association["varList"->varList],
    Association["var"->var]
  ]
]


Options[readVAR]={"ltrim"->True,"imagic"->{}};
readVAR[sim_,iVAR_,OptionsPattern[]]:=Module[{nprocx,nprocy,nprocz,varList,varAll,var},
  {nprocx,nprocy,nprocz}=nProc[sim];
  varList=varName[sim]//Keys;
  
  varAll=ConstantArray[{},{nprocx,nprocy,nprocz}];
  Monitor[
    Do[
      var=readVARProc[sim,iproc,iVAR,{{nprocx,nprocy,nprocz},varList},
        "ltrimInner"->True,"ltrimOuter"->OptionValue["ltrim"],
        "imagic"->OptionValue["imagic"]
      ];
      varAll[[Sequence@@(1+var["ipxyz"])]]=var["var"],
      {iproc,Range[0,nprocx*nprocy*nprocz-1]}
    ],StringJoin["Reading proc: ",ToString@(1+iproc)," / ",ToString[nprocx*nprocy*nprocz]]
  ];
  
  (* varAll is of size (nprocx,nprocy,nprocz,nvar,sx,sy,sz) *)
  (* where sxyz can be mxyz-ghxyz or mxyz-2ghxyz depending on ipxyz *)
  varAll=Flatten[varAll,{{4},{1,5},{2,6},{3,7}}];
  
  Join[
  Association["t"->var["t"]],
  AssociationThread[var["varList"]->varAll]
  ]
]


(* ::Section:: *)
(*Interpolation*)


pcIntpltVar[var_,vars_List]:=Table[
  Interpolation@Transpose[{Transpose[var/@{"x","y","z"}],var[f]}],
{f,vars}]


(* ::Chapter:: *)
(*End*)


End[]


Protect[
  readVARNProc,readVARN,
  readVARProc,readVAR,
  tSnap,
  pcIntpltVar
]


EndPackage[]
