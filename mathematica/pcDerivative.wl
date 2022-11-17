(* ::Package:: *)

(* :Name: pcDerivative` *)
(* :Author: Hongzhe Zhou, Stockholm, 2021*)
(* :Version: 0.1 *)
(* :Summary:
    Contains 6th-order finite difference functions.
*)


BeginPackage["pcDerivative`","pcReadBasic`"]


(* ::Chapter:: *)
(*Usage messages*)


der::usage="der[f,mx,my,mz,ghi,di,i,n] returns the nth derivative w.r.t. i of f.
Input:
  f:  1D List of length mx*my*mz. Must be sorted by coordiantes
  ghi: number of ghost cells
  di:  grid spacing
  i:  1, 2, or 3 as x, y, or z derivative
Output:
  \!\(\*SuperscriptBox[\(d\), \(n\)]\)f/\!\(\*SuperscriptBox[\(di\), \(n\)]\) as a mx*my*mz length List"

der2::usage="der2[f,mx,my,mz,{ghi,ghj},{di,dj},{i,j}] returns d^2f/didj.
Input:
  f:  1D List of length mx*my*mz. Must be sorted by coordiantes
  ghi,ghij: numbers of ghost cells
  di,dj:  grid spacing
  i,j:  1, 2, or 3 as x, y, or z derivative
Output:
  d^2f/didj as a mx*my*mz length List"

grad::usage="grad[f,mx,my,mz,ghx,ghy,ghz,dx,dy,dz] returns the gradient of f.
Input:
  f:  1D List of length mx*my*mz. Must be sorted by coordiantes
  ghx,ghy,ghz: numbers of ghost cells
  dx,dy,dz:  grid spacings
Output:
  {df/dx,df/dy,df/dz}, each as an mx*my*mz length List"

div::usage="div[{fx,fy,fz},mx,my,mz,ghx,ghy,ghz,dx,dy,dz] returns the divergence of f.
Alternatively, div[sim,var,vec] directly takes divergence on a var file.
Input:
  fx,fy,fz: Three components of f, each as a 1D List of length mx*my*mz.
             Must be sorted by coordiantes.
  ghx,ghy,ghz: Numbers of ghost cells.
  dx,dy,dz: Grid spacings.
  sim: Directory of the run.
  var: The return of readVARN[...]. Must include ghost zones.
  vec: The three components of the vector, e.g., {\"uu1\",\"uu2\",\"uu3\"}.
Output:
  div(f), as an mx*my*mz length List"

curl::usage="curl[{fx,fy,fz},mx,my,mz,ghx,ghy,ghz,dx,dy,dz] returns the curl of f.
Input:
  fx,fy,fz: three components of f, each as a 1D List of length mx*my*mz.
             Must be sorted by coordiantes
  ghx,ghy,ghz: numbers of ghost cells
  dx,dy,dz:  grid spacings
Output:
  curl(f), as three mx*my*mz length Lists"

curl2::usage="curl2[{fx,fy,fz},mx,my,mz,ghx,ghy,ghz,dx,dy,dz] returns the curl of curl of f.
Input:
   fx,fy,fz: three components of f, each as a 1D List of length mx*my*mz.
             Must be sorted by coordiantes
  ghx,ghy,ghz: numbers of ghost cells
  dx,dy,dz:  grid spacings
Output:
  curl(curl(f)), as three mx*my*mz length Lists"

laplacian::usage="laplacian[f,mx,my,mz,ghx,ghy,ghz,dx,dy,dz] returns the laplacian of f.
Input:
  f:  1D List of length mx*my*mz. Must be sorted by coordiantes
  ghx,ghy,ghz: numbers of ghost cells
  dx,dy,dz:  grid spacings
Output:
  laplacian(f), as three mx*my*mz length Lists"

gradDiv::usage="gradDiv[{fx,fy,fz},mx,my,mz,ghx,ghy,ghz,dx,dy,dz] returns grad(div f).
Input:
  fx,fy,fz: three components of f, each as a 1D List of length mx*my*mz.
             Must be sorted by coordiantes
  ghx,ghy,ghz: numbers of ghost cells
  dx,dy,dz:  grid spacings
Output:
  grad(div f), as three mx*my*mz length Lists"

trim::usage="trim[sim,f] trims the ghost zones of data f.
Input:
  sim: Directory of the run.
  f: A List of length mx*my*mz."


Begin["`Private`"]


(* ::Chapter:: *)
(*Cartesian coordinates*)


(* ::Section:: *)
(*nth derivatives (6th order)*)


(*1st derivatives*)
dd1D[mi_,gh_,1]:=With[{a=3/4,b=-3/20,c=1/60},
  SparseArray[{
    Band[{1,2}]->a,Band[{2,1}]->-a,
    Band[{1,3}]->b,Band[{3,1}]->-b,
    Band[{1,4}]->c,Band[{4,1}]->-c,
    Band[{1,-1-2gh},{1,-1-2gh}]->-a,
    Band[{1,-2-2gh},{2,-1-2gh}]->-b,
    Band[{1,-3-2gh},{3,-1-2gh}]->-c,
    Band[{-1,1+2gh}]->a,
    Band[{-2,1+2gh}]->b,
    Band[{-3,1+2gh}]->c
  },{mi,mi}]
]

(*2nd derivatives*)
dd1D[mi_,gh_,2]:=With[{a=3/2,b=-3/20,c=1/90,d=-49/18},
  SparseArray[{
    Band[{1,1}]->d,
    Band[{1,2}]->a,Band[{2,1}]->a,
    Band[{1,3}]->b,Band[{3,1}]->b,
    Band[{1,4}]->c,Band[{4,1}]->c,
    Band[{1,-1-2gh},{1,-1-2gh}]->a,
    Band[{1,-2-2gh},{2,-1-2gh}]->b,
    Band[{1,-3-2gh},{3,-1-2gh}]->c,
    Band[{-1,1+2gh}]->a,
    Band[{-2,1+2gh}]->b,
    Band[{-3,1+2gh}]->c
  },{mi,mi}]
]

ddx[f_,mx_,my_,mz_,ghx_,dx_,n_Integer]/;n>=1:=With[{ddx1D=dd1D[mx,ghx,n]},
  1/dx^n*Flatten[ddx1D . Partition[f,my*mz]]
]

ddy[f_,mx_,my_,mz_,ghy_,dy_,n_Integer]/;n>=1:=With[{ddy1D=dd1D[my,ghy,n]},
  Module[{tmp},
    tmp=ArrayReshape[Transpose[ArrayReshape[f,{mx,my,mz}],1<->2],{my,mx*mz}];
    1/dy^n*Flatten[Transpose[ArrayReshape[ddy1D . tmp,{my,mx,mz}],1<->2]]
  ]
]

ddz[f_,mx_,my_,mz_,ghz_,dz_,n_Integer]/;n>=1:=With[{ddz1D=dd1D[mz,ghz,n]},
  Module[{tmp},
    tmp=ArrayReshape[Transpose[ArrayReshape[f,{mx,my,mz}],1<->3],{mz,mx*my}];
    1/dz^n*Flatten[Transpose[ArrayReshape[ddz1D . tmp,{mz,my,mx}],1<->3]]
  ]
]

(*wrapper for nth derivatives w.r.t. the same variable*)
der[f_,mx_,my_,mz_,ghi_,di_,i_,n_Integer]/;n>=1:=Switch[i,
  1,ddx[f,mx,my,mz,ghi,di,n],
  2,ddy[f,mx,my,mz,ghi,di,n],
  3,ddz[f,mx,my,mz,ghi,di,n]
]

(*wrapper for 2nd derivatives*)
der2[f_,mx_,my_,mz_,{ghi_,ghi_},{di_,di_},{i_,i_}]:=der[f,mx,my,mz,ghi,di,i,2]
der2[f_,mx_,my_,mz_,{ghi_,ghj_},{di_,dj_},{i_,j_}]/;i!=j:=Module[{tmp},
  tmp=der[f,mx,my,mz,ghi,di,i,1];
  der[tmp,mx,my,mz,ghj,dj,j,1]
]


(* ::Section:: *)
(*Derivative operators*)


grad[f_,mx_,my_,mz_,ghx_,ghy_,ghz_,dx_,dy_,dz_]:=MapThread[
  der[f,mx,my,mz,#1,#2,#3,1]&,{{ghx,ghy,ghz},{dx,dy,dz},{1,2,3}}
]

div[{fx_,fy_,fz_},mx_,my_,mz_,ghx_,ghy_,ghz_,dx_,dy_,dz_]:=MapThread[
  der[#1,mx,my,mz,#2,#3,#4,1]&,{{fx,fy,fz},{ghx,ghy,ghz},{dx,dy,dz},{1,2,3}}
]//Total

curl[{fx_,fy_,fz_},mx_,my_,mz_,ghx_,ghy_,ghz_,dx_,dy_,dz_]:={
  der[fz,mx,my,mz,ghy,dy,2,1]-der[fy,mx,my,mz,ghz,dz,3,1],
  der[fx,mx,my,mz,ghz,dz,3,1]-der[fz,mx,my,mz,ghx,dx,1,1],
  der[fy,mx,my,mz,ghx,dx,1,1]-der[fx,mx,my,mz,ghy,dy,2,1]
}

laplacian[f_,mx_,my_,mz_,ghx_,ghy_,ghz_,dx_,dy_,dz_]:=MapThread[
  der[f,mx,my,mz,#1,#2,#3,2]&,{{ghx,ghy,ghz},{dx,dy,dz},{1,2,3}}
]//Total

gradDiv[{fx_,fy_,fz_},mx_,my_,mz_,ghx_,ghy_,ghz_,dx_,dy_,dz_]:=MapThread[Plus[
      der2[fx,mx,my,mz,{ghx,#1},{dx,#2},{1,#3}],
      der2[fy,mx,my,mz,{ghy,#1},{dy,#2},{2,#3}],
      der2[fz,mx,my,mz,{ghz,#1},{dz,#2},{3,#3}]
    ]&,{{ghx,ghy,ghz},{dx,dy,dz},{1,2,3}}]

curl2[{fx_,fy_,fz_},mx_,my_,mz_,ghx_,ghy_,ghz_,dx_,dy_,dz_]:=Module[{lapf},
  lapf=laplacian[#,mx,my,mz,ghx,ghy,ghz,dx,dy,dz]&/@{fx,fy,fz};
  gradDiv[{fx,fy,fz},mx,my,mz,ghx,ghy,ghz,dx,dy,dz]-lapf
]


(* ::Section:: *)
(*Derivatives on VAR files*)


trim[sim_,ff_]:=Module[{mx,my,mz,gh1,gh2,gh3,tmp},
  {mx,my,mz,gh1,gh2,gh3}=readDim[sim]/@{"mx","my","mz","gh1","gh2","gh3"};
  tmp=ff//ArrayReshape[#,{mx,my,mz}]&;
  tmp[[gh1+1;;mx-gh1,gh2+1;;my-gh2,gh3+1;;mz-gh3]]//Flatten
]

div[sim_,var_,vec_List]:=Module[{mx,my,mz,gh1,gh2,gh3,dxyz,tmp},
  {mx,my,mz,gh1,gh2,gh3}=readDim[sim]/@{"mx","my","mz","gh1","gh2","gh3"};
  dxyz=var/@{"dx","dy","dz"};
  trim[sim,div[var/@vec,mx,my,mz,gh1,gh2,gh3,Sequence@@dxyz]]
]

grad[sim_,var_,sc_]:=Module[{mx,my,mz,gh1,gh2,gh3,dxyz,tmp},
  {mx,my,mz,gh1,gh2,gh3}=readDim[sim]/@{"mx","my","mz","gh1","gh2","gh3"};
  dxyz=var/@{"dx","dy","dz"};
  trim[sim,#]&/@grad[var[sc],mx,my,mz,gh1,gh2,gh3,Sequence@@dxyz]
]


(* ::Chapter:: *)
(*End*)


End[]


Protect[
  der,der2,
  grad,div,curl,curl2,laplacian,
  trim
]


EndPackage[]
