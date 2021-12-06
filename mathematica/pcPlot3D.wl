(* ::Package:: *)

(* :Name: pcPlot3D` *)
(* :Author: Hongzhe Zhou, Stockholm, 2021*)
(* :Version: 0.1 *)
(* :Summary:
    Functions that make 3D plots.
*)


BeginPackage["pcPlot3D`","pcReadBasic`","pcRead2D`","pcPlot`"]


(* ::Chapter:: *)
(*Usage messages*)


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
  An Image object."

makeBoxes::usage="makeBoxes[T,L,R,minmax:,plotStyle:] reads video files and returns a list of
  3D box view of var. You can then make the legend by hand using, e.g., DensityPlot[...].
  Currently it is assumed that top at zmax, left at ymin, and right at xmax.
  If not, you can change the ticks by hand by saying Ticks->{...}.
  For non-cubic boxes, simply say, for example, BoxRatios->{1,1,2}, and correspondingly change AxesLabel; then
  you will get a tall box, etc.
Input:
  T,L,R: Lists. Must be the returns of readSlices[...].
  minmax: List of length 2. Specifies the plotting range. If not present, then the MinMax of the
          data from all the three surfaces will be used.
  plotStyle: Rule, Optional. Will overwrite other options of the plot.
Output:
  A list of Image objects."


Begin["`Private`"]


(* ::Chapter:: *)
(*Functions*)


(* ::Section:: *)
(*Make a 3D box snapshot*)


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
      Axes->True,AxesLabel->{
        \!\(\*
TagBox[
StyleBox["\"\<\\!\\(\\*StyleBox[\\\"x\\\",FontSlant->\\\"Italic\\\"]\\)/\\[Pi]\>\"",
ShowSpecialCharacters->False,
ShowStringCharacters->True,
NumberMarks->True],
FullForm]\),
        \!\(\*
TagBox[
StyleBox["\"\<\\!\\(\\*StyleBox[\\\"y\\\",FontSlant->\\\"Italic\\\"]\\)/\\[Pi]\>\"",
ShowSpecialCharacters->False,
ShowStringCharacters->True,
NumberMarks->True],
FullForm]\),
        \!\(\*
TagBox[
StyleBox["\"\<\\!\\(\\*StyleBox[\\\"z\\\",FontSlant->\\\"Italic\\\"]\\)/\\[Pi]\>\"",
ShowSpecialCharacters->False,
ShowStringCharacters->True,
NumberMarks->True],
FullForm]\)
      },
      Ticks->ConstantArray[{{0,-1},{0.5,0},{1,1}},3],
      LabelStyle->pcLabelStyle,Frame->True,FrameStyle->pcLabelStyle,
      ViewPoint->{2,-2,1.7},ImageSize->Medium,
      Epilog->{Inset[Style[
          "\!\(\*StyleBox[\"t\",\nFontSlant->\"Italic\"]\)="<>ToString@time,pcLabelStyle],
        Scaled[{0.9,0.94}]]}
    ]//Rasterize[#,ImageSize->500,RasterSize->1000]&
  ]

makeBoxes[slT_,slL_,slR_,minmax_List:{},plotStyle___Rule]:=Module[{top,left,right,t,minMax,i},
  {top,t}=Most@slT;
  {left,right}=First/@{slL,slR};
  minMax=If[minmax!={},minmax,MinMax[Flatten[{top,left,right}]]];
  Print["PlotRange is {All,All,All,",minMax,"}"];
  SetSharedVariable[i];
  i=0;
  Monitor[
    ParallelMap[(i++;makeBox[Sequence@@#,minmax,plotStyle])&,Transpose[{top,left,right,t}]],
    StringJoin["Making box ",ToString@i,"/",top//Length//ToString]
  ]
]


(* ::Chapter:: *)
(*End*)


End[]


Protect[
  makeBox,makeBoxes
]


EndPackage[]
