(* ::Package:: *)

(* :Name: pcPowerSpec` *)
(* :Author: Hongzhe Zhou, Stockholm, 2022*)
(* :Version: 0.1 *)
(* :Summary:
    Some functions to compute power spectra.
*)


BeginPackage["pcPowerSpec`"]


(* ::Chapter:: *)
(*Usage messages*)


power1D::usage="power1D[data,{nx,ny,nz}] computes the 1D power spectrum of data. The normalization is
Total[powerSpectrum]=Mean[data^2].
Input:
  data: 1D List of dimension {nx,ny,nz}: { f111, f112, f113, ..., f121, f122, ... }.
Output:
  The power spectrum as a List of length Max[nx,ny,nz]/2, corresponding to wavenumbers 0,1,2,..."

power1DintZ::usage="power1DintZ[data,nx,Options] first does an integration over z, then
computes the 1D power spectrum in the xy plane. The normalization is Total[powerSpectrum]=Total[data^2].
Input:
  data: 1D List of length nx^3. By default it assumes the data corresponds to the coordinates
        {{x1,y1,z1},{x2,y1,z1},...}, but this can be changed by the option \"lxyz\".
Output:
  The power spectrum as a List of length nx/2, corresponding to wavenumbers 0,1,2,...,nx/2-1.
Options:
  \"lxyz\": By default \"xyz\". This can be changed to be \"zyx\" to instruct that the data
            is in the order {{x1,y1,z1},{x1,y1,z2},...}."


Begin["`Private`"]


(* ::Chapter:: *)
(*Functions*)


(* ::Section:: *)
(*1D power spectrum*)


(* shell-integrated power spectrum *)
power1D[data_,{nx_,ny_,nz_}]:=Module[{data1,fft,nmax,kx,kxyz,knorm,kshell,pos,power,fact},
  (* put data into a cube, do fft, and back to a 1D List *)
  data1=ArrayReshape[data,{nx,ny,nz}];
  fft=Flatten@Fourier[data1];
  
  (* kx, ky, and kz coordinates *)
  nmax=Max[nx,ny,nz];
  kx=N@Range[0,nmax-1]/.x_/;x>(nmax/2-1):>x-nmax;
  kxyz=Outer[List,kx,kx,kx]//Flatten[#,2]&;
  knorm=Norm/@kxyz;
  
  (* kr, and determine which positions belongs to each kr *)
  kshell=Range[0,nmax/2-1];
  pos=Position[knorm,x_/;Round[x]==#]&/@kshell;
  
  (* shell integration and normalization*)
  power=Re[Total[Extract[fft*Conjugate[fft],#]]&/@pos];
  Return[power/Total[power]*Mean[data^2]]
]


(* first integrate over z, then compute 1d power spectrum in the xy plane *)
Options[power1DintZ]={"lxyz"->"xyz"};
power1DintZ[data_,nx_,OptionsPattern[]]:=Module[{xy,data1,fft,kx,kxy,knorm,kshell,pos,power},
  xy=Switch[OptionValue["lxyz"],
    (* {{x1,y1},{x2,y1},...} *)
    "xyz",Table[{i,j},{k,nx},{j,nx},{i,nx}]//Flatten[#,2]&,
    (* {{x1,y1},{x1,y2},...} *)
    "zyx",Table[{i,j},{i,nx},{j,nx},{k,nx}]//Flatten[#,2]&
  ];
  (* Gather data by their xy coordinates *)
  data1={xy,data}//Transpose//GatherBy[#,First]&;
  (* Integration over z; now the data is alwyas in the {{x1,y1},{x1,y2},...} order *)
  data1=Last/@Sort[Total/@data1];
  (* Put data back into a squre, do fft, and back to a 1D List again *)
  data1=Partition[data1,nx];
  fft=Flatten@Fourier[data1];
  
  (* kx and ky coordiantes *)
  kx=N@Range[0,nx-1]/.x_/;x>(nx/2-1):>x-nx;
  kxy=Outer[List,kx,kx]//Flatten[#,1]&;
  knorm=Norm/@kxy;
  (* kr, and determine which positions belongs to each kr *)
  kshell=Range[0,nx/2-1];
  pos=Position[knorm,x_/;Round[x]==#]&/@kshell;

  (* shell integration and normalization*)
  power=Re[Total[Extract[fft*Conjugate[fft],#]]&/@pos];
  Return[power/Total[power]*Total[data^2]]
]


(* ::Chapter:: *)
(*End*)


End[]


Protect[
  power1D,power1DintZ
]


EndPackage[]
