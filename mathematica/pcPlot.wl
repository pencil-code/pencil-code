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

makeBox::usage="makeBox[dataTop,dataLeft,dataRight,time,minmax:{},plotStyle:] visualizes data
  in a 3D box.
Input:
  dataTop: List. The values to be plotted on the top plane. Must have depth=3.
  dataLeft: List. The values to be plotted on the left plane. Must have depth=3.
  dataRight: List. The values to be plotted on the right plane. Must have depth=3.
  time: NumericQ. Time will be displayed on the top right in the figure.
  minmax: List, Optional. Must have the form {min,max}, which will then determine
          the color range. If not provided then MinMax[{dataTop,dataLeft,dataRight}//Flatten]
          will be used.
  plotStyle: Rule, Optional. Will overwrite other options of the plot.
Output:
  A Graphics3D object."

makeBoxes::usage="makeBoxes[sim,var,plotStyle:] reads video files and returns a list of
  3D box view of var. You can then make the legend by hand using, e.g., DensityPlot[...].
Input:
  sim: String. The run directory.
  var: String. The variable needed to plot. Must be present in video.in.
  plotStyle: Rule, Optional. Will overwrite other options of the plot.
Output:
  A list of Graphics3D object."


Begin["`Private`"]


(* ::Chapter:: *)
(*Functions*)


(* ::Section:: *)
(*Plot style related*)


pcLabelStyle=Directive[Thick,Black,14,FontFamily->"Times"];
pcPlotStyle[]:={
  Map[SetOptions[#,{
      PlotRange->All,Frame->True,LabelStyle->pcLabelStyle,
      FrameStyle->pcLabelStyle,ImageSize->{300,300/GoldenRatio}
    }]&,{
      Plot,LogPlot,LogLogPlot,LogLinearPlot,DensityPlot,
      ListPlot,ListLogPlot,ListLogLogPlot,ListLogLinearPlot,ListLinePlot,
      ListDensityPlot
    }
  ];
  Map[SetOptions[#,{
      PlotLegends->Automatic,ColorFunction->"Rainbow"
    }]&,{DensityPlot,ListDensityPlot}
  ];
  Map[SetOptions[#,{
      InterpolationOrder->0
    }]&,{ListDensityPlot}
  ];
}


(* ::Section:: *)
(*Generate a 2D video from video files*)


Options[showVideo]={"lFrames"->False};
showVideo[sim_,var_,loc_,plotStyle_List:{},aniStyle_List:{},OptionsPattern[]]:=
  Module[{sp,slice,time,pos,minmax,frames},
    sp=", "<>Which[
      StringMatchQ[loc,"xy"~~___],"z=",
      StringMatchQ[loc,"yz"~~___],"x=",
      StringMatchQ[loc,"xz"~~___],"y=",
      True,""
    ];
  {slice,time,pos}=readSlice[sim,var,loc];
  minmax=MinMax[slice//Flatten];
  frames=MapThread[ListDensityPlot[#1,plotStyle,
    ColorFunction->"Rainbow",ColorFunctionScaling->False,
    PlotLegends->BarLegend[{"Rainbow",minmax}],FrameTicks->None,
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


makeBox[dataTop_,dataLeft_,dataRight_,time_,minmax_List:{},plotStyle___Rule]:=
  Module[{minMax,cf,mkImg,img,mkSlice},
    minMax=If[minmax!={},minmax,MinMax[Flatten[{dataTop,dataLeft,dataRight}]]];
    cf=ColorData["Rainbow"][Rescale[#,minMax]]&;
    mkImg[data_]:=ListDensityPlot[data,PlotRange->All,
      PlotLegends->None,Frame->False,ColorFunction->cf,ColorFunctionScaling->False,
      ImagePadding->None,PlotRangePadding->None,ImageMargins->None,ImageSize->{200,200}
    ];
    {img["T"],img["L"],img["R"]}=mkImg/@{dataTop,dataLeft,dataRight};
    mkSlice[sp_]:=Module[{pol},
      pol=Polygon[Switch[sp,
          "T",{{1,0,1},{1,1,1},{0,1,1},{0,0,1}},
          "L",{{1,0,0},{1,0,1},{0,0,1},{0,0,0}},
          "R",{{1,1,0},{1,1,1},{1,0,1},{1,0,0}}
        ],VertexTextureCoordinates->{{1,0}, {1,1}, {0,1}, {0,0}}];
      Graphics3D[{Texture[img[sp]],pol}]
    ];
    Show[mkSlice["T"],mkSlice["L"],mkSlice["R"],plotStyle,
      PlotRange->{{0,1},{0,1},{0,1}},Lighting->{"Ambient",White},
      AxesLabel->{"x/(2\[Pi])","y/(2\[Pi])","z/(2\[Pi])"},Axes->True,
      LabelStyle->pcLabelStyle,Frame->True,FrameStyle->pcLabelStyle,
      ViewPoint->{2,-2,1.7},ImageSize->Medium,
      Epilog->{Inset[Style["t="<>ToString@time,pcLabelStyle],Scaled[{0.9,0.94}]]}
    ]
  ]

makeBoxes[sim_,var_,plotStyle___Rule]:=Module[{top,left,right,t,minmax,legend},
  {top,t}=Most@readSlice[sim,var,"xy"];
  {left,right}=First[readSlice[sim,var,#]]&/@{"xz","yz"};
  minmax=MinMax[Flatten[{top,left,right}]];
  Print["PlotRange is {All,All,All,",minmax,"}"];
  MapThread[
    makeBox[#1,#2,#3,#4,minmax,plotStyle,PlotLabel->var]&,
    {top,left,right,t}
  ]
]


(* ::Chapter:: *)
(*End*)


End[]


Protect[
  pcLabelStyle,pcPlotStyle,
  showVideo,
  showSlice,showSliceVector,
  makeBox,makeBoxes
]


EndPackage[]
