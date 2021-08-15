(* ::Package:: *)

(* :Name: pcParticleStalker` *)
(* :Author: Hongzhe Zhou, Stockholm, 2021*)
(* :Version: 0.1 *)
(* :Summary:
    Functions for stalker particles.
*)


BeginPackage["pcParticleStalker`","pcReadBasic`"]


(* ::Chapter:: *)
(*Usage messages*)


nStalkerSnap::usage="nStalkerSnap[sim] returns the number of snapshots of stalker particles."

stalkerHead::usage="stalkerHead[sim] returns a list of quantities carried by stalker particles."

nStalker::usage="nStalker[sim] returns the number of stalker particles."

readStalker::usage="readStalker[sim] reads all stalker particle snapshots on all profiles.
Input:
  sim: String. Directory of the run.
Output:
  { {t1,t2,...},\[IndentingNewLine]    { {xp of particles at t1},{xp of particles at t2},... },\[IndentingNewLine]    { {yp of particles at t1},{yp of particles at t2},... },\[IndentingNewLine]    ...\[IndentingNewLine]  }"


Begin["`Private`"]


(* ::Chapter:: *)
(*Functions*)


(* ::Section:: *)
(*Read basic information*)


nStalkerSnap[sim_]:=Import[sim<>"/data/tstalk.dat"][[1,2]]

stalkerHead[sim_]:=StringSplit[Import[sim<>"/data/particles_stalker_header.dat"][[1,1]],","]

nStalker[sim_]:=Import[sim<>"/data/pdim.dat"][[1,3]]


(* ::Section:: *)
(*Read one round of particles_stalker.dat from one proc*)


readStalkerProc1[file_,pre_,nhead_]:=Module[{start,tsp,nv,id,sdata},
  start=BinaryRead[file,"Integer32"];
  tsp=BinaryRead[file,pre];
  nv=BinaryRead[file,"Integer32"];
  BinaryRead[file,"Integer32"];
  
  If[nv<1, {"t"->tsp}//Association//Return];

  BinaryRead[file,"Integer32"];
  id=BinaryRead[file,ConstantArray["Integer32",nv]];
  BinaryRead[file,"Integer32"];

  BinaryRead[file,"Integer32"];
  sdata=BinaryRead[file,ConstantArray[pre,nv*nhead]];
  BinaryRead[file,"Integer32"];
  
  {"t"->tsp,Rule@@@Transpose[{id,Partition[sdata,nhead]}]}//Flatten//Association
]


(* ::Section:: *)
(*Read particles_stalker.dat from one proc*)


readStalkerProc[sim_,proc_]:=With[
  {file=sim<>"/data/proc"<>ToString[proc]<>"/particles_stalker.dat",
  pre=readDim[sim]["precision"],
  head=stalkerHead[sim]},
  Module[{nhead=Length[head],nsnap=nStalkerSnap[sim]},
    Close[file]//Quiet;
    
    (*returns length nsnap Association, each points to an Association id\[Rule]pvars*)
    Table[i->readStalkerProc1[file,pre,Length[head]],{i,nsnap}]//Association
  ]
]
(*
  Example: readStalkerProc[sim,0][2] gives the 2nd snapshot on proc0
  <|"t"\[Rule]t,id1\[Rule]{...},...|>
*)


(* ::Section:: *)
(*Read particles_stalker.dat from all procs*)


readStalker[sim_]:=With[
  {nproc=nProc[sim],
  nsnap=nStalkerSnap[sim],
  nsp=nStalker[sim],
  header=Thread[#->Range[Length@#]]&@stalkerHead[sim]//Association},
  Module[{data},
    (*data[iproc][isnap]=id\[Rule]vars*)
    data=Table[i->readStalkerProc[sim,i],{i,0,nproc-1}]//Association;

    Table[With[{tmp=Join@@Table[data[iproc][isnap],{iproc,0,nproc-1}]},
      {tmp["t"],Sequence@@Transpose[tmp/@Range[nsp]]}
    ],{isnap,nsnap}]//Transpose
  ]
]


(* ::Chapter:: *)
(*End*)


End[]


Protect[
  nStalkerSnap,stalkerHead,nStalker,
  readStalker
]


EndPackage[]
