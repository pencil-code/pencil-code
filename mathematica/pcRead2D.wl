(* ::Package:: *)

(* :Name: pcRead2D` *)
(* :Author: Hongzhe Zhou, Stockholm, 2021*)
(* :Version: 0.1 *)
(* :Summary:
    Reads 2D slice files.
*)


BeginPackage["pcRead2D`","pcReadBasic`"]


(* ::Chapter:: *)
(*Usage messages*)


readStride::usage="readStride[sim,slice] reads stride files.
Input:
  sim: String. Directory of the run.
  slice: String, which slice to read. E.g., xy, xy2.
Output:
  An Association object of the min/max/step indices for each directions.
"
readSlice::usage="readSlice[sim,var,slice] reads slice files.
Input:
  sim: String. Directory of the run
  var: String. Variable to be read
  slice: String, which slice to read
Output:
  {slices,times,Union[positions]}
Example:
  readSlice[sim,\"uu1\",\"xy2\"]"

read2DAve::usage="read2DAve[sim,sp] reads 1-dimensionally averaged 2-dimensional
files, i.e., those variables specified in {xyz}aver.in.
Input:
  sim: String. Directory of the run.
  sp: String, which direction is averaged out: x, y, or z.
Output:
  An Association object of time, grid, and fields."


Begin["`Private`"]


(* ::Chapter:: *)
(*Functions*)


(* ::Section:: *)
(*Read slice files (video.in)*)


readStride[sim_,sl_]:=Module[{dim,dir,n,file,stride},
  dim=readDim[sim];
  dir=StringTake[sl,#]&/@{{1},{2}};
  n=dim["n"<>#]&/@dir;
  
  file=FileNameJoin[{sim,"data/stride_"<>StringTake[sl,2]<>".dat"}];
  If[!FileExistsQ[file],
    (* default values if stride file does not exist *)
    stride={1,n[[1]],1,1,n[[2]],1},
    
    (* read stride file *)
    stride=Flatten@Import[file];
    Switch[Length[stride],
      2, stride={1,n[[1]],stride[[1]],1,n[[2]],stride[[2]]},
      6, stride=stride,
      _, Print["stride file bad structure. Failed."];Return@$Failed]
  ];
  
  stride=stride/.{"nxgrid"->dim["nx"],"nygrid"->dim["ny"],"nzgrid"->dim["nz"]};
  
  AssociationThread[
    Flatten[Table[i<>j,{i,dir},{j,{"min","max","step"}}]]->stride
  ]
]

readSlice[sim_,var_,sl_]:=Module[
  {file,p,stride,n1,n2,readOneTime,data,slices,times,positions},
  file=sim<>"/data/slice_"<>var<>"."<>sl;
  readSlice::nofile=StringJoin[file," does not exist. Please check."];
  readSlice::badformat=StringJoin[file," has bad format. Please check."];
  If[!FileExistsQ[file],Message[readSlice::nofile];Return[$Failed]];
  
  (* determine data dimension *)
  p=readDim[sim]["precision"];
  stride=readStride[sim,sl]//Values;
  n1=Length[Range@@stride[[1;;3]]];
  n2=Length[Range@@stride[[4;;6]]];
  
  (* read *)
  readOneTime[x_]:=Module[{nstart=Last[x],slice,t,pos,nend,nnext},
    slice=BinaryRead[file,ConstantArray[p,n1*n2]];
    t=BinaryRead[file,p];
    pos=BinaryRead[file,p];
    nend=BinaryRead[file,"Integer32"];
    If[nstart!=nend,Message[readSlice::badformat];Return[$Failed]];
    Return[{Partition[slice,n1],t,pos,BinaryRead[file,"Integer32"]}]
  ];
  Close[file]//Quiet;
  data=NestWhileList[readOneTime,{BinaryRead[file,"Integer32"]},(Last[#]=!=EndOfFile)&];
  {slices,times,positions}=Most[Transpose[Rest@data]];
  {slices,times,Union[positions]}
]


(* ::Section:: *)
(*1D-averaged 2D files ({xyz}aver.in)*)


read2DAve[sim_,sp_]:=Module[{
    sp1,sp2,i,j,
    nsnap,procs,vars,
    file,n1,n2,pre,grid,t,f,out
  },
  
  (* error and warning messages *)
  read2DAve::noend="Warning: File did not reach EndofFile for proc `1`";
  read2DAve::difft="Error: inconsistent times on processors. A list of proc id vs. time is returned.";
  
  (* determine which coordinates on the plane *)
  {sp1,sp2}=DeleteCases[{"x","y","z"},sp];
  {i,j}=Flatten[Position[{"x","y","z"},#]&/@{sp1,sp2}];
  
  (* global variables to all procs *)
  (* number of snapshots *)
  nsnap=Import[sim<>"/data/t2davg.dat"][[1,-1]]-1;
  (* processors which have the data *)
  procs=Cases[Range[0,Times@@(nProc[sim])-1],x_/;readDim[sim,x]["ip"<>sp]==0];
  (* variables *)
  vars=Cases[Flatten@Import[sim<>"/"<>sp<>"aver.in"],x_/;StringTake[x,1]!="#"];
  
  (* read data from each processor *)
  Do[
    file=sim<>"/data/proc"<>ToString[iproc]<>"/"<>sp<>"averages.dat";
    Close[file]//Quiet;
    {n1,n2,pre}=readDim[sim,iproc]/@{"n"<>sp1,"n"<>sp2,"precision"};
    (* grid on this processor *)
    grid=Outer[List,Sequence@@(readGrid[sim,iproc][[{i,j}]]),1]//Transpose//Flatten[#,1]&;
    
    (* loop over chunks of snapshots *)
    t[iproc]={};
    f[iproc]=Module[{tmp},Table[
        BinaryRead[file,"Integer32"];
        t[iproc]={t[iproc],BinaryRead[file,pre]};  (* time *)
        BinaryRead[file,"Integer32"];
        
        BinaryRead[file,"Integer32"];
        tmp=BinaryRead[file,ConstantArray[pre,n1*n2*Length[vars]]];  (* fields *)
        BinaryRead[file,"Integer32"];
        
        Join[{grid},Partition[tmp,n1*n2]]
        
        ,nsnap]
      ]; (* end of reading from this processor *)
      t[iproc]=Flatten[t[iproc]];
      
      (* consistency check *)
      If[BinaryRead[file,"Integer32"]=!=EndOfFile,Message[read2DAve::noend,iproc]];
      Close[file]
  ,{iproc,procs}];
  
  (* check if time on each processor agree *)
  If[Not[Equal@@(t/@procs)],
    Message[read2DAve::difft];Return[Transpose[{procs,t/@procs}]]
  ];
  
  (* combine data into an Association *)
  out=Apply[Join,Transpose[f/@procs,{4,2,1,3}],{2}];
  Join[
    AssociationThread[vars->Rest@out],
    Association["grid"->Transpose[out[[1,1]]]],
    Association["t"->t[procs[[1]]]]
  ]
]


(* ::Chapter:: *)
(*End*)


End[]


Protect[
  readStride,readSlice,read2DAve
]


EndPackage[]
