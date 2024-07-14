(* ::Package:: *)

(* :Name: pcPlot3D` *)
(* :Author: Hongzhe Zhou, Stockholm, 2021*)
(* :Version: 0.1 *)
(* :Summary:
    Functions that make 3D plots.
*)


BeginPackage["pcPlot3D`","pcPlot`"]


(* ::Chapter:: *)
(*Usage messages*)


mk3DBox::usage="mk3DBox[dataLeft,dataRight,dataTop] visualizes data in a 3D box.
Input:
  dataLeft, dataRight, dataTop: The data to be plotted on the left, right, and top
                                planes, respectively. These shall be 2D Lists.
Options:
  \"PlotStyle\": The options to be inherited by ArrayPlot for making the three images.
  \"BoxStyle\": The options to be inherited by Show for making the 3D box.
Example:
  {l,r,t}=readSlice[sim,\"lnrho\",#][[1,-1]]&/@{\"xz\",\"yz\",\"xy\"};
  mk3DBox[l,r,t,\"PlotStyle\"->{ColorFunction->pcColors[\"Rainbow\"]},\"BoxStyle\"->{Background->Black}]"


Begin["`Private`"]


(* ::Chapter:: *)
(*Functions*)


(* ::Section:: *)
(*Make a 3D box snapshot*)


Options[mk3DBox]={"PlotStyle"->{ImageSize->Automatic},"BoxStyle"->{ImageSize->Automatic}};
mk3DBox[dataLeft_,dataRight_,dataTop_,OptionsPattern[]]:=Module[{minmax,img,img3D},
  minmax={-1,1}*Max[{dataLeft,dataRight,dataTop}//Flatten//Abs];
  
  {img["L"],img["R"],img["T"]}=ArrayPlot[
    Reverse@Transpose[#],OptionValue["PlotStyle"],
    ColorFunction->pcColors["BlueGreenYellow"],
    ImagePadding->None,
    Frame->None,
    PlotRange->minmax
  ]&/@{dataLeft,dataRight,dataTop};
  
  With[{tmp1={{0,1},{1,1},{1,0},{0,0}},tmp2={{0,0}, {1,0}, {1,1}, {0,1}}},
    img3D["L"]= Graphics3D[{Texture[img["L"]],
      Polygon[{#1,0,#2}&@@@tmp1,
      VertexTextureCoordinates->tmp2]
    }];
    img3D["R"]= Graphics3D[{Texture[img["R"]],
      Polygon[{1,#1,#2}&@@@tmp1,
      VertexTextureCoordinates->tmp2]
    }];
    img3D["T"]= Graphics3D[{Texture[img["T"]],
      Polygon[{#1,#2,1}&@@@tmp1,
      VertexTextureCoordinates->tmp2]
    }];
  ];
  
  Show[img3D["L"],img3D["R"],img3D["T"],
    OptionValue["BoxStyle"],
    PlotRange->{{0,1},{0,1},{0,1}},
    Lighting->{"Ambient",White}, ViewPoint->{2,-2,1.7},
    Frame->None,Axes->None
  ]
]


(* ::Chapter:: *)
(*End*)


End[]


Protect[
  mk3DBox
]


EndPackage[]
