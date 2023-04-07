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


pcStreamLineSolve2D::usage="pcStreamLineSolve2D[intpFunc,pt0,n,options] calculates a field line in 2D.
Inputs:
  intpFunc: the InterpolatingFunction of a 2D vector field, i.e., intpFunc[x,y]={fx,fy}. Should be
            consistent with the \"Coordinates\" option, i.e., if \"Coordinates\"->\"Spherical\", then
            x and y are the radial and colatitude coordinates, respectively, and intpFunc should also
            produce the radial and colatitute components of the vector field.
  pt0:      List of length 2. The starting point of the field line. Should be consistent with the
            \"Coordinates\" option.
  n:        Real. The length of the field line. If negative, then back-trace the field line.
Options:
  \"Coordinates\": Specifies the coordinate system. Possible values are: \"Spherical\".
                   \"Cartesian\" is coded but not tested.
  \"BoundCond\":   Specifies the boundary conditions. Possible values are:
                   \"Spherical\": the field line stops at any of the boundaries.
                   \"pp\":        periodic in both x and y. Coded but not tested.
                   \"sp\":        shear-periodic in x and periodic in y. Coded but not tested.
  \"deltay\":      The shift in y when using the shear-periodic boundary condition.
Output:
  A List of length 2 consisting of two InterpolatingFunctions. These are the Cartesian (regardless of
  \"Coordinates\") positions of the field line point. Each of them takes in one argument which is the
  affine parameter. Typically its largest (smallest, if n<0) possible value is smaller (larger, if n<0)
  than n. This can be checked by saying streamLine[[1]][\"Domain\"] if streamLine is the output of this
  function."

pcStreamLineData2D::usage="pcStreamLineData2D[intpFunc,pt0,n,nData:128, options] calculates a field line
in 2D that is ready to plot.
Inputs:
  intpFunc, pt0, n: The same as those in pcStreamLineSolve2D.
  nData: Integer, Optional. Specifies how many data points to use along the field line.
Options:
  The same as those in pcStreamLineSolve2D.
Output:
  A List of length <=nData."

pcStreamLineSolve::usage="This is intended for plotting 3D field lines, but outdated and should be updated."


Begin["`Private`"]


(* ::Chapter:: *)
(*Functions*)


Options[pcStreamLineSolve2D]={"Coordinates"->"Cartesian","BoundCond"->"pp","deltay"->0};
pcStreamLineSolve2D[intpFunc_,pt0_List,n_,OptionsPattern[]]:=Module[{
  coordBounds1,coordBounds2,coordToXy,xyToCoord,fieldToXy,
  testBC, applyBC, testAndApplyBC,
  fx,fy,dxyds,eqns,x,y,s
  },
  
  (* set up coordinate system *)
  coordBounds1=intpFunc["Domain"][[1]];
  coordBounds2=intpFunc["Domain"][[2]];
  Switch[OptionValue["Coordinates"],
    "Cartesian",
      coordToXy[coord1_,coord2_]:={coord1,coord2};
      xyToCoord[x_,y_]:={x,y};
      fieldToXy[{coord1_,coord2_},{f1_,f2_}]:={f1,f2},
    "Spherical",
      (* coords = radial and colatitude *)
      coordToXy[coord1_,coord2_]:=coord1*{Sin[coord2],Cos[coord2]};
      xyToCoord[x_,y_]:={Norm[{x,y}],ArcCos[y/Norm[{x,y}]]};
      fieldToXy[{coord1_,coord2_},{f1_,f2_}]:={f1*Sin[coord2]+f2*Cos[coord2],f1*Cos[coord2]-f2*Sin[coord2]},
    _,
      Print["Coordinates="<>OptionValue["Coordinates"]<>" not implemented!"];Return[$Failed]
  ];
  
  (* periodic boundary? *)
  Switch[OptionValue["BoundCond"],
    "pp",
      testBC[coord1_,coord2_]:=True;
      applyBC[coord1_,coord2_]:=N@MapThread[#1-Floor[#1-#2,#3-#2]&,
        {{coord1,coord2},Transpose[{coordBounds1,coordBounds2}]//Sequence}],
    "sp",
      testBC[coord1_,coord2_]:=True;
      applyBC[coord1_,coord2_]:=N@MapThread[#1-Floor[#1-#2,#3-#2]&,
        {{coord1,coord2+OptionValue["deltay"]},Transpose[{coordBounds1,coordBounds2}]//Sequence}],
    "Spherical",
      testBC[coord1_,coord2_]:=Between[coord1,coordBounds1] && Between[coord2,coordBounds2];
      applyBC[x_,y_]:=N@{x,y},
    _,
      Print["BoundCond="<>OptionValue["BoundCond"]<>" not implemented!"];Return[$Failed]
  ];
  If[Not[testBC@@pt0],Print["The given initial point at ",pt0," is outside the boundaries!"];Return[$Failed]];
  
  (* Take in (x,y), transform to (coord1,coord2). Then test and/or apply b.c., and transform back to xy *)
  testAndApplyBC[x_,y_]:=With[{coords=xyToCoord[x,y]},
    If[testBC@@coords,coordToXy@@(applyBC@@coords),"outside"]
  ];
  
  (* field values in Cartesian coordinates *)
  fx["outside"]:=0.;
  fy["outside"]:=0.;
  fx[{x_?NumericQ, y_?NumericQ}]:=With[{coords=xyToCoord[x,y]},fieldToXy[coords,intpFunc@@coords][[1]]];
  fy[{x_?NumericQ, y_?NumericQ}]:=With[{coords=xyToCoord[x,y]},fieldToXy[coords,intpFunc@@coords][[2]]];
  
  (* we solve the equation (x'(s),y'(s))=(fx(xyz),fx(xyz))/Sqrt[fx^2+fx^2] *)
  (* also allow for n<0, i.e., back-tracing *)
  dxyds[x_,y_]:=With[{xy=testAndApplyBC[x,y]},Normalize@Through[{fx,fy}[xy]]];
  eqns=Join[
    Thread[{x'[s],y'[s]}==Sign[n]*dxyds[x[s],y[s]]],
    Thread[{x[0],y[0]}==coordToXy@@pt0]
  ];
  Head/@NDSolveValue[eqns,{x[s],y[s]},{s,0,Abs[n]}]
]

pcStreamLineData2D[intpFunc_,pt0s_List,n_,ndata_Integer:128,opts:OptionsPattern[]]:=Module[{stlines,domain,data},
  stlines=pcStreamLineSolve2D[intpFunc,#,n,opts]&/@pt0s;
  data=Table[If[Head/@l==={InterpolatingFunction,InterpolatingFunction},
      domain=l[[1]]["Domain"]//Flatten; Through[l[#]]&/@Subdivide[Sequence@@domain,ndata],
      $Failed
    ],
    {l,stlines}
  ];
  Print["Fail to compute the following lines:"];
  Print[Extract[pt0s,Position[data,$Failed]]];
  DeleteCases[data,$Failed]
]

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
  pcStreamLineSolve2D,pcStreamLineData2D,
  pcStreamLineSolve
]


EndPackage[]
