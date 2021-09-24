(* ::Package:: *)

(* :Name: pcPlot` *)
(* :Author: Hongzhe Zhou, Stockholm, 2021*)
(* :Version: 0.1 *)
(* :Summary:
    Some plot style.
*)


BeginPackage["pcPlot`","pcReadBasic`","pcRead2D`"]


(* ::Chapter:: *)
(*Usage messages*)


pcLabelStyle::usage="Font and size.";
pcPlotStyle::usage="Set some plot styles.";

showVideo::usage="showVideo[sim,var,loc,plotStyle:{},aniStyle:{},\"lFrames\"->False] generates a video from video files.
Input:
  sim: String. The run directory
  var: String. The variable to be read. Must be present in video.in
  loc: String. Which slice to read; e.g., \"xy2\" or \"yz\"
  plotStyle: List. Optional. Additional options for ListDensityPlot
  aniStyle: List. Optional. Additional options for ListAnimate
  \"lFrames\": Option. False by default. If ->True then returns a list of frames.
Output:
  An animation."

showSlice::usage="showSlice[data,var,{sp,loc},plotStyle:{},\"ldata\"->False] plots a slice.
Input:
  data: Must be from readVARN[sim,iVAR]
  var: String. Variable to read; e.g. \"uu1\"
  sp: String. \"x\", \"y\", or \"z\"
  loc: Integer, the number of the slice
  plotStyle: Optional. Changes plotting style of ListDensityPlot
  \"ldata\": Option. False by default. If ->True then returns data on the slice.
Output:
  A density plot of var at location loc in the sp direction"

showSliceVector::usage="showSliceVector[data,var,{sp,loc},plotStyle:{}] plots vectors
  within a slice; e.g, for a z=0 slice only the x and y components are plotted.
Input:
  data: Must be from readVARN[sim,iVAR]
  var:  String. Variable to read; e.g. \"uu\" or \"bbb\"
  sp:  String. \"x\", \"y\", or \"z\"
  loc:  Integer, the number of the slice
  plotStyle: Optional. Changes plotting style of ListDensityPlot
Output:
  A vector plot of var at location loc in the sp direction"

(*some problems with using Hyperplane*)
(*renderBox::usage="renderBox[sim,var,xy,yz,xz,{time},{style}] gives 3D plots.
Input:
  sim: String. Directory of the run
  var: String, variable to be read
  xy,yz,xz: String. Optional, top, left and right slice specifications;
            by default xy=\"xy2\", yz=\"yz\", xz=\"xz\"
  time: List. one or more values. The plots are made at time closest to the specified ones.
  style: List. Optional, additional style specification
Output:
  A List of 3D plots"*)


Begin["`Private`"]


(* ::Chapter:: *)
(*Functions*)


(* ::Section:: *)
(*Plot style related*)


pcColorFunction[z_]:=Blend[{RGBColor[0.368417, 0.506779, 0.709798],White,RGBColor[0.880722, 0.611041, 0.142051]},z]
pcColorFunction[x_,y_,z_]:=pcColorFunction[z]

pcLabelStyle=Directive[Thick,Black,14,FontFamily->"Times"];
pcPlotStyle[]:={
  Map[SetOptions[#,{
      PlotRange->All,Frame->True,LabelStyle->pcLabelStyle,
      FrameStyle->pcLabelStyle,ImageSize->{300,300/GoldenRatio}
    }]&,{
      Plot,LogPlot,LogLogPlot,LogLinearPlot,
      ListPlot,ListLogPlot,ListLogLogPlot,ListLogLinearPlot,ListLinePlot
    }
  ];
  Map[SetOptions[#,{
      PlotRange->All,Frame->True,LabelStyle->pcLabelStyle,PlotLegends->Automatic,
      FrameStyle->pcLabelStyle,ImageSize->200,ColorFunction->pcColorFunction
    }]&,{ListDensityPlot}
  ];
}


(* ::Section:: *)
(*Generate a 2D video from video files*)


Options[showVideo]={"lFrames"->False};
showVideo[sim_,var_,loc_,plotStyle_List:{},aniStyle_List:{},OptionsPattern[]]:=
  Module[{sp,slice,time,pos,max,cf,frames},
    sp=", "<>Which[
      StringMatchQ[loc,"xy"~~___],"z=",
      StringMatchQ[loc,"yz"~~___],"x=",
      StringMatchQ[loc,"xz"~~___],"y=",
      True,""
    ];
  {slice,time,pos}=readSlice[sim,var,loc];
  max=1.1*Max[slice//Flatten//Abs];
  cf[z_]:=pcColorFunction[(z+max)/(2*max)];
  frames=MapThread[ListDensityPlot[#1,plotStyle,
    ColorFunction->cf,ColorFunctionScaling->False,
    PlotLegends->BarLegend[{cf,{-max,max}}],FrameTicks->None,
    PlotLabel->StringJoin["t=",ToString[NumberForm[#2,3]],sp,ToString[NumberForm[pos[[1]],3]]]
  ]&,{slice,time}];
  If[OptionValue["lFrames"],frames,ListAnimate[frames,aniStyle,AnimationRunning->False]]
]


(* ::Section:: *)
(*Plot a slice from a readVARN[sim,iVAR] data*)


Options[showSlice]={"ldata"->False};
showSlice[data_,var_String,{sp_String,loc_Integer},plotStyle_List:{},OptionsPattern[]]:=
  Module[{f,x1,x2,pos,r,lx,ly,lz},
    {lx,ly,lz}=data/@{"lx","ly","lz"};
    f=data[var];
    Switch[sp,
      "x",{x1,x2}=data/@{"y","z"};pos={loc,;;,;;};r=lz/ly,
      "y",{x1,x2}=data/@{"x","z"};pos={;;,loc,;;};r=lz/lx,
      "z",{x1,x2}=data/@{"x","y"};pos={;;,;;,loc};r=ly/lx
    ];
    {f,x1,x2}=Flatten[
      Part[ArrayReshape[#,{lx,ly,lz}],Sequence@@pos]
    ]&/@{f,x1,x2};
    If[OptionValue["ldata"],Transpose[{x1,x2,f}]//Return];
    ListDensityPlot[Transpose[{x1,x2,f}],
      Join[plotStyle,{
        PlotLegends->Automatic,FrameStyle->Directive[Black,Thick],
        AspectRatio->r,
        LabelStyle->Directive[Thick,Black,14,FontFamily->"Times"]
      }]
    ]
  ]

Options[showSliceVector]={"data"->False};
showSliceVector[data_,var_String,{sp_String,loc_Integer},plotStyle_List:{},OptionsPattern[]]:=
  Module[{f1,f2,x1,x2,ff,xx,pos,r,lx,ly,lz},
    {lx,ly,lz}=data/@{"lx","ly","lz"};
    Switch[sp,
      "x",{f1,f2,x1,x2}=data/@{var<>"2",var<>"3","y","z"};pos={loc,;;,;;};r=lz/ly,
      "y",{f1,f2,x1,x2}=data/@{var<>"1",var<>"3","x","z"};pos={;;,loc,;;};r=lz/lx,
      "z",{f1,f2,x1,x2}=data/@{var<>"1",var<>"2","x","y"};pos={;;,;;,loc};r=ly/lx
    ];
    {f1,f2,x1,x2}=Flatten[
      Part[ArrayReshape[#,{lx,ly,lz}],Sequence@@pos]
    ]&/@{f1,f2,x1,x2};
    ff=Transpose[{f1,f2}];
    xx=Transpose[{x1,x2}];
    If[OptionValue["data"],Transpose[{xx,ff}]//Return];
    ListVectorPlot[Transpose[{xx,ff}],
      Join[plotStyle,{
        FrameStyle->Directive[Black,Thick], AspectRatio->r,
        LabelStyle->Directive[Thick,Black,14,FontFamily->"Times"]
      }]
    ]
  ]


(* ::Section:: *)
(*3D Plots*)


renderBox[sim_,var_,xy_String:"xy2",yz_String:"yz",xz_String:"xz",{time__},style_List:{}]:=
  Module[{nx,ny,nz,slices,l,left,right,top,t,pos,plot},
    {nx,ny,nz}=readDim[sim]/@{"nx","ny","nz"};
    t=readSlice[sim,var,xy][[2]];
    Print["Available snapshots at: ",t];
    pos=Nearest[t->"Index",#][[1]]&/@{time};
    slices=readSlice[sim,var,#]&/@{xy,yz,xz};
    Print["{x,y,z} position: ",slices[[#,3]]&/@{2,3,1}];
    Do[
      {top,left,right}=(#[[1,i]])&/@slices;
      l=ConstantArray[0.,{nz,ny,nx}];
      l[[nz,All,All]]=top;
      l[[All,1,All]]=right;
      l[[All,All,1]]=left;
      plot[i]=ListSliceDensityPlot3D[l,
        {Hyperplane[{0,0,1},{1,1,nz}],Hyperplane[{0,1,0},{1,1,1}],Hyperplane[{1,0,0},{1,1,1}]},
        Flatten[{style,
         BoxRatios->{nx,ny,nz},Axes->None,Boxed->False,
         ViewPoint->{-0.7nx,-0.6ny,0.32nz},ViewVertical->{0,0,1},
         ImageSize->{400,400},PlotLegends->Automatic,
         PlotLabel->"\!\(\*
StyleBox[\"t\",\nFontSlant->\"Italic\"]\)="<>ToString[t[[i]]],
         LabelStyle->Directive[Thick,Black,14,FontFamily->"Times"]
       }]
      ],{i,pos}
    ];
    plot/@pos
  ]


(* ::Chapter:: *)
(*End*)


End[]


Protect[
  pcLabelStyle,pcPlotStyle,
  showVideo,
  showSlice,showSliceVector
]


EndPackage[]
