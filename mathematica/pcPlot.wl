(* ::Package:: *)

(* :Name: pcPlot` *)
(* :Author: Hongzhe Zhou, Stockholm, 2021*)
(* :Version: 0.1 *)
(* :Summary:
    Define general plot styles, and functions that make 2D plots.
*)


BeginPackage["pcPlot`","pcReadBasic`","pcRead1D`","pcRead2D`"]


(* ::Chapter:: *)
(*Usage messages*)


pcLabelStyle::usage="Font and size.";
pcPlotStyle::usage="Set some plot styles.";
pcPopup::usage="Using DisplayFunction->pcPopup will make a plot in a pop-up window."

spaceTimeDiag::usage="spaceTimeDiag[sim,plane,var,plotStyle:] makes a butterfly diagram from planar averaged data.
Input:
  sim: String. Directory of the simulation.
  plane: String. \"xy\" or \"yz\" or \"xz\", specifying which plane is averaged over.
  var: String. Must be present in \"xyaver.in\" etc.
  plotStyle: Optional. Can be several Rule objects. Will overwrite the default plot styles.
Output:
  A space-time diagram with 1/3 aspect ratio, and time coordinates in code units. The space cooridnates
are normalized to (-1,1). Also print out the MinMax values of the data for the possibility of making
legends by hand."

showSlice::usage="showSlice[data,var,{sp,loc},plotStyle:] plots a slice.
Input:
  data: Must be from readVARN[sim,iVAR].
  var: String. Variable to read; e.g. \"uu1\".
  sp: String. \"x\", \"y\", or \"z\".
  loc: Numeric. Will plot the slice closest to this position.
  plotStyle: Optional. Will overwrite the current style of ListDensityPlot.
Output:
  A density plot of var at location loc in the sp direction"

showSliceVector::usage="showSliceVector[data,var,{sp,loc},plotStyle:] plots a vector field
  in the slice spcified by sp and loc. The in-plane componets are shown as arrows, and the
  out-of-plane componet is shown in a density map.
Input:
  data: Must be from readVARN[sim,iVAR]
  var: String. Variable to read; e.g. \"uu\" or \"bbb\"
  sp: String. \"x\", \"y\", or \"z\"
  loc: Numeric. Will plot the slice closest to this position.
  plotStyle: Optional. Will overwrite the current plotting style.
Output:
  A vector plot of var at location loc in the sp direction"

pcPanelDataPlot::usage="pcPanelDataPlot[data,{m,n},style:{},Options] makes a plot of m times n panels,
with shared axes.
Input:
  data: List. Each element should be a data set like {{x1,y1},{x2,y2},...}.
  {m,n}: m columns and n rows. Better to have m*n=Length[data].
  style: List. Optional. User-specified plot styles for all panels.
Options:
  \"ImageSize0\" by default 300. The ImageSize for each individual panel.
  \"ImagePadding0\" by default {45,6}. The ImagePadding for each individual panel.
      The first number specifies the ones with frame labels, and the second number
      specifies the ones without frame labels.
  \"Spacing\" by default {0.2,0.5}. The horizontal and vertical spacings between panels."

pcPanelPlot::usage="pcPanelPlot[plots,{m,n},Options] makes a plot of m times n panels,
with shared axes.
Input:
  plots: A List of plots. Can already have FrameTicks.
  {m,n}: m columns and n rows. Better to have m*n=Length[plots].
Options:
  \"ImageSize0\" by default 300. The ImageSize for each individual panel.
  \"ImagePadding0\" by default 25. The ImagePadding for each individual panel.
  \"Spacing\" by default {0.2,0.5}. The horizontal and vertical spacings between panels."


Begin["`Private`"]


(* ::Chapter:: *)
(*Functions*)


(* ::Section:: *)
(*Plot style related*)


pcLabelStyle=Directive[AbsoluteThickness[1],Black,14,FontFamily->"Times"];
pcPlotStyle[]:=Module[{setOps},
  setOps[ops_List,funcs_List]:=Map[SetOptions[#,ops]&,funcs];
  (*General options for all plots*)
  setOps[{
      PlotRange->All,Frame->True,LabelStyle->pcLabelStyle,
      FrameStyle->pcLabelStyle,ImageSize->{300,300/GoldenRatio},
      ImagePadding->{{60,10},{60,10}}
    },{
      Plot,LogPlot,LogLogPlot,LogLinearPlot,DensityPlot,
      ListPlot,ListLogPlot,ListLogLogPlot,ListLogLinearPlot,ListLinePlot,
      ListDensityPlot,ListVectorPlot
    }];
  (*Options for 1D ListPlot's*)
  setOps[{
      Joined->True
    },{
      ListPlot,ListLogPlot,ListLogLogPlot,ListLogLinearPlot,ListLinePlot
    }];
  (*Options for 2D plots*)
  setOps[{
      PlotLegends->Automatic,ColorFunction->"Rainbow"
    },{
      DensityPlot,ListDensityPlot
    }];
  (*Options for ListDensity Plot*)
  setOps[{
      InterpolationOrder->0
    },{
      ListDensityPlot
    }];
]

pcPopup[plot_]:=CreateDocument[plot,
  "CellInsertionPointCell"->Cell,ShowCellBracket->False,
  WindowElements->{},WindowFrame->"Generic",WindowSize->All,
  WindowTitle->None,WindowToolbars->{}
]


(* ::Section:: *)
(*Butterfly diagram*)


spaceTimeDiag[sim_,sl_,var_,plotStyle___Rule]:=Module[{t,f,nt,nx,gf},
  PrintTemporary["Reading data..."];
  {t,f}=readAves[sim,sl,var];
  {nt,nx}=Dimensions[f];

  gf=f[[1;;-1;;Ceiling[nt/nx/3],;;]]//Transpose;
  Print["The MinMax of data is ",gf//Flatten//MinMax];

  PrintTemporary["Making plots..."];
  ListDensityPlot[gf,plotStyle,
    DataRange->{t[[{1,-1}]],{-1,1}},
    AspectRatio->1/3,InterpolationOrder->0
  ]
]


(* ::Section:: *)
(*Plot a slice from a readVARN[sim,iVAR] data*)


showSlice[data_,var_String,{sp_String,loc_?NumericQ},plotStyle___Rule]:=
  Module[{f,x1,x2,x3,pos,plotData},
    f=data[var];
    {x1,x2,x3}=data/@Switch[sp,
        "x",{"y","z","x"},"y",{"x","z","y"},"z",{"x","y","z"}
      ];
    pos=Nearest[x3->"Index",loc];
    Print["Plotting the slice at ",sp,"=",x3[[pos//First]]];
    plotData=Transpose[Extract[#,List/@pos]&/@{x1,x2,f}];
    ListDensityPlot[plotData,plotStyle,
      AspectRatio->Abs[(Subtract@@MinMax[x2])/(Subtract@@MinMax[x1])]
    ]
  ]

(*plot from VAR data*)
showSliceVector[data_Association,var_String,{sp_String,loc_?NumericQ},plotStyle___Rule]:=
  Module[{f1,f2,f3,x1,x2,x3,r,pos,x12,v12,vecData,denData},
    {f1,f2,f3}=data[var<>ToString@#]&/@Switch[sp,
      "x",{2,3,1},"y",{1,3,2},"z",{1,2,3}
    ];
    {x1,x2,x3}=data/@Switch[sp,
      "x",{"y","z","x"},"y",{"x","z","y"},"z",{"x","y","z"}
    ];
    r=Abs[(Subtract@@MinMax[x2])/(Subtract@@MinMax[x1])];
    pos=Nearest[x3->"Index",loc];
    Print["Plotting the slice at ",sp,"=",x3[[pos//First]]];
    Print["Plotting range of the out-of-plane component: ",MinMax[f3]];
    x12=Transpose[Extract[#,List/@pos]&/@{x1,x2}];
    v12=Transpose[Extract[#,List/@pos]&/@{f1,f2}];
    vecData=Transpose[{x12,v12}];
    denData=Flatten/@Transpose[{x12,Extract[f3,List/@pos]}];
    Show[
      ListDensityPlot[denData],
      ListVectorPlot[vecData,VectorScaling->"Linear",
        VectorColorFunction->None,VectorStyle->White],
      plotStyle,AspectRatio->r,ImageSize->Medium
    ]
]


(* ::Section:: *)
(*Plot multiple panels with common axes*)


Options[pcPanelDataPlot]={"ImageSize0"->300,"ImagePadding0"->25,"Spacingx"->-2}
pcPanelDataPlot[data_,{m_,n_},style_List:{},OptionsPattern[]]:=Module[{
  nData,lpos,shift,xRange,yRange,styleAll,
  tkStyle,imgPd,imgSz,imgSz0,imgPd0,spx},
  pcPanelDataPlot::nonrect="Warning: The bottom ticks may show incorrectly.";
  pcPanelDataPlot::insuffpanel="Error: Insufficient number of panels specified.";
  nData=Length[data];
  If[m*n<nData,Message[pcPanelDataPlot::insuffpanel];Return[$Failed]];
  If[m*n!=nData,Message[pcPanelDataPlot::nonrect]];
  
  imgSz0=OptionValue["ImageSize0"];
  imgPd0=OptionValue["ImagePadding0"];
  spx=OptionValue["Spacingx"];
  lpos[i_]:=List[
    (*left and right*)
    {MemberQ[Range[1,m*n,m],i],MemberQ[Range[m,m*n,m],i]},
    (*bottom and top*)
    {MemberQ[Range[m*n-m+1,m*n],i],(*MemberQ[Range[1,m],#]*)False}
  ];
  
  (*overall plot styles*)
  shift[{x1_,x2_}]:={x1-0.1(x2-x1),x2+0.1(x2-x1)};
  xRange=data[[;;,;;,1]]//Flatten//MinMax//shift;
  yRange=data[[;;,;;,2]]//Flatten//MinMax//shift;
  styleAll={style,ImageSize->imgSz0,PlotRange->{xRange,yRange},FrameTicks->True,FrameTicksStyle->pcLabelStyle};
  
  (*individual plot styles*)
  tkStyle[i_]:=(FrameTicksStyle->Map[Directive[FontOpacity->#]&,lpos[i]/.{True->1,False->0},{2}]);
  imgPd[i_]:=(ImagePadding->(lpos[i]/.{True->imgPd0,False->1}));
  imgSz[i_]:=Module[{w=imgSz0,h=imgSz0/GoldenRatio},
    If[!lpos[i][[1,1]],w=w-imgPd0];If[!lpos[i][[1,2]],w=w-imgPd0];
    If[!lpos[i][[2,1]],h=h-imgPd0];If[!lpos[i][[2,2]],h=h-imgPd0];
    ImageSize->{w,h}
  ];
  
  Table[ListPlot[data[[i]],tkStyle[i],imgPd[i],imgSz[i],styleAll],{i,nData}]//Partition[#,UpTo@m]&//Grid[#,Spacings->{spx,0.5}]&
]


Options[pcPanelPlot]={"ImageSize0"->300,"ImagePadding0"->{45,6},"Spacing"->{0.2,0.5}}
pcPanelPlot[plots_,{m_,n_},OptionsPattern[]]:=Module[{
  nData,lpos,ladd,tkStyle,imgPd,imgSz,imgSz0,imgPd0,spxy},
  pcPanelPlot::insuffpanel="Error: Insufficient number of panels specified.";
  nData=Length[plots];
  If[m*n<nData,Message[pcPanelPlot::insuffpanel];Return[$Failed]];
  ladd=If[m*n!=nData,True,False];
  
  imgSz0=OptionValue["ImageSize0"];
  imgPd0=OptionValue["ImagePadding0"];
  spxy=OptionValue["Spacing"];
  lpos[i_]:=List[
    (*left and right*)
    {MemberQ[Range[1,m*n,m],i],MemberQ[Range[m,m*n,m],i]},
    (*bottom and top*)
    {MemberQ[Range[m*n-m+1,m*n],i]||And[ladd,MemberQ[Range[nData-m+1,m*(n-1)],i]],
    (*MemberQ[Range[1,m],#]*)False}
  ];
  
  (*individual plot styles*)
  tkStyle[i_]:=(FrameTicksStyle->Map[Directive[FontOpacity->#]&,lpos[i]/.{True->1,False->0},{2}]);
  imgPd[i_]:=(ImagePadding->(lpos[i]/.{True->imgPd0[[1]],False->imgPd0[[2]]}));
  imgSz[i_]:=Module[{w=imgSz0,h=imgSz0/GoldenRatio},
    If[!lpos[i][[1,1]],w=w-imgPd0];If[!lpos[i][[1,2]],w=w-imgPd0];
    If[!lpos[i][[2,1]],h=h-imgPd0];If[!lpos[i][[2,2]],h=h-imgPd0];
    ImageSize->{w,h}
  ];
  
  Table[Show[plots[[i]],tkStyle[i],imgPd[i],imgSz[i]],{i,nData}]//Partition[#,UpTo@m]&//Grid[#,Spacings->spxy,Alignment->Top]&
]


(* ::Chapter:: *)
(*End*)


End[]


Protect[
  pcLabelStyle,pcPlotStyle,pcPopup,
  spaceTimeDiag,
  showVideo,
  showSlice,showSliceVector,
  pcPanelDataPlot,pcPanelPlot
]


EndPackage[]
