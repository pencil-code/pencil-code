(* ::Package:: *)

(* :Name: pcPlot` *)
(* :Author: Hongzhe Zhou, Shanghai, 2022*)
(* :Version: 0.1 *)
(* :Summary:
    This module takes care of plotting stream lines of vector fields.
*)


BeginPackage["pcStreamLine`","pcReadVAR`"]


(* ::Chapter:: *)
(*Usage messages*)


pcStreamLineSolve::usage=""


Begin["`Private`"]


(* ::Chapter:: *)
(*Functions*)


Options[pcStreamLineSolve]={"Boundcond"->"ppp","xyz0"->{-\[Pi],-\[Pi],-\[Pi]},"xyz1"->{\[Pi],\[Pi],\[Pi]},"deltay"->0};
pcStreamLineSolve[vars_List,pt0_List,n_,OptionsPattern[]]:=Module[{applyBC,xyz0,xyz1,dxyzds,eqns,x,y,z,s},
  xyz0=OptionValue["xyz0"];
  xyz1=OptionValue["xyz1"];
  
  Switch[OptionValue["Boundcond"],
    "ppp", applyBC[x_,y_,z_]:=N@MapThread[#1-Floor[#1-#2,#3-#2]&,{{x,y,z},xyz0,xyz1}],
    "spp", applyBC[x_,y_,z_]:=N@MapThread[#1-Floor[#1-#2,#3-#2]&,{{x,y+OptionValue["deltay"],z},xyz0,xyz1}],
    _,     Print["Boundary condition ot implemented."];Return[$Failed]
  ];
  
  (* we solve the equation (x'(s),y'(s),z'(s))=(f1(xyz),f2(xyz),f3(xyz))/Sqrt[f1^2+f2^2+f3^2] *)
  dxyzds[x_,y_,z_]:=With[{xyz2=applyBC[x,y,z]},Normalize@Through[vars@@xyz2]];  
  eqns=Join[
    Thread[{x'[s],y'[s],z'[s]}==dxyzds[x[s],y[s],z[s]]],
    Thread[{x[0],y[0],z[0]}==pt0]
  ];
  Head/@NDSolveValue[eqns,{x[s],y[s],z[s]},{s,0,n}]
]


(* ::Chapter:: *)
(*End*)


End[]


Protect[
  pcStreamLineSolve
]


EndPackage[]
