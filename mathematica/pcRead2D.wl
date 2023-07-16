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


Begin["`Private`"]


(* ::Chapter:: *)
(*Functions*)


(* ::Section:: *)
(*Read slice files*)


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
      6, ,
      _, Print["stride file bad structure. Failed."];Return@$Failed]
    ];
    
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


(* ::Chapter:: *)
(*End*)


End[]


Protect[
  readStride,readSlice
]


EndPackage[]
