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


readSlice[sim_,var_,sl_]:=
  With[{
      file=sim<>"/data/slice_"<>var<>"."<>sl,
      dvid=readParamNml[sim,"run.in","DVID"],
      tlast=readTS[sim,"t"][[-1]]
    },
    readSlice::nofile=StringJoin[file," does not exist. Please check."];
    readSlice::badformat=StringJoin[file," has bad format. Please check."];
    Module[{p,nx,ny,nz,n1,n2,readOneTime,data,slices,times,positions},
      If[!FileExistsQ[file],Message[readSlice::nofile];Return[$Failed]];
      {p,nx,ny,nz}=readDim[sim]/@{"precision","nx","ny","nz"};
      {n1,n2}=Switch[StringTake[sl,2],"xy",{nx,ny},"yz",{ny,nz},"xz",{nx,nz}];
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
  ]


(* ::Chapter:: *)
(*End*)


End[]


Protect[
  readSlice
]


EndPackage[]
