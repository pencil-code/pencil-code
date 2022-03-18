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
      ListDensityPlot,ListVectorPlot
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


(* ::Chapter:: *)
(*End*)


End[]


Protect[
  pcLabelStyle,pcPlotStyle,pcPopup,
  spaceTimeDiag,
  showVideo,
  showSlice,showSliceVector
]


EndPackage[]
