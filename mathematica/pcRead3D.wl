(* ::Package:: *)

(* :Name: pcRead3D` *)
(* :Author: Hongzhe Zhou, Stockholm, 2021*)
(* :Version: 0.1 *)
(* :Summary:
    Reads 3D data.
*)


BeginPackage["pcRead3D`","pcReadBasic`","pcRead1D`"]


(* ::Chapter:: *)
(*Usage messages*)


read3DFFT::usage="read3DFFT[sim,sp] reads 3D data from file data/fft*.dat. The dimension of the data
 should be given by kout_max in run.in.
Input:
  sim: Run directory
  sp: uu or oo
Output:
  { {t1,t2,...},
    {{k1,k2,k3},...},
    {List of data}
  }"


Begin["`Private`"]


(* ::Chapter:: *)
(*Functions*)


read3DFFT[sim_,file_]:=Module[{km,kvec,t,f},
  km=readParamNml[sim,"run.in","KOUT_MAX"];
  kvec=Flatten[Table[{kx,ky,kz},{kz,-km,km},{ky,-km,km},{kx,-km,km}],2];
  {t,f}=read1D[sim,file];
  {t,kvec,f}
]


(* ::Chapter:: *)
(*End*)


End[]


Protect[
  read3DFFT
]


EndPackage[]
