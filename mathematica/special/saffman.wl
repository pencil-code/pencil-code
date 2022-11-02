(* ::Package:: *)

(* :Name: saffman` *)
(* :Summary:
    Defines getParam[] which produces various useful parameters.
*)


BeginPackage["saffman`"]
Unprotect["saffman`*"];
ClearAll["saffman`*"];
ClearAll["saffman`Private`*"];


(* ::Chapter:: *)
(*Usage messages*)


pDecay::usage="Computes the decaying exponent t^(-p), as a function of t.";
lB::usage="Computes scale of B."
pqDiag::usage="pq diagram."


Begin["`Private`"]


(* ::Chapter:: *)
(*Functions*)


Needs["pc`"]


pDecay[t_,x_,i_:1]:=-i*Most[t/x]*Differences[x]/Differences[t]

lB[sim_String,f_:Identity]:=Module[{t,spec,k},
  {t,spec}=read1D[sim,"power_mag.dat"]//Transpose//f//Transpose;
  spec=Rest/@spec;
  k=Range[1,getParam[sim,"Nx"]/2-1];
  Map[Total[#/k]/Total[#]&,spec]
]
lB[sim_String,t_List,spec_List]:=Module[{k},
  k=Range[1,getParam[sim,"Nx"]/2-1];
  Map[Total[#/k]/Total[#]&,Rest/@spec]
]


Options[pqDiag]={"lnotPQ"->False,"lPanelLabel"->""};
pqDiag[sim_,p_,q_,t_,style_List:{},OptionsPattern[]]:=Module[{pq,style2,n,r},
  style2={AxesOrigin->{0,0},ImagePadding->{{50,10},{50,10}},ImageSize->360};
  pq=ListPlot[List/@(Transpose[{q,p}]),
    Joined->False,
    PlotStyle->Red,
    (*PlotMarkers->(Graphics[{Disk[]},ImageSize->8(1+#)]&/@Rescale[t])*)
    PlotMarkers->Table[{"\[EmptyCircle]",s},{s,20/1.3*(0.3+#)&/@Rescale[t]}]
  ];
  If[OptionValue["lnotPQ"],
    Show[pq,style2,style],
    n=2getParam[sim,"nHyperEta"];
    r=readParamNml[sim,"run.in","ETA_TDEP_EXPONENT"];
    Show[pq,
      Plot[2(1-x),{x,0,1},PlotStyle->Directive[AbsoluteThickness[1],Black]],
      Plot[(2n+2r)/(n-1)-(4n-2)/(n-1)*x,{x,0,1},PlotStyle->Directive[Black,Dashed]],
      ListPlot[Table[{x,2.5 x},{x,0,1,0.5}],PlotStyle->Directive[Black,Dotted]],
      ListPlot[Table[Table[{x,(1+k) x},{x,0,1,0.5}],{k,Range[0,4]}],PlotStyle->Directive[Black,Dotted]],
      PlotRange->{{0,1},{0,2}},
      Epilog->{
        Inset[Style["\[Beta]=4",pcLabelStyle],Scaled[{0.3,0.92}]],
        Inset[Style["\[Beta]=3",pcLabelStyle],Scaled[{0.55,0.92}]],
        Inset[Style["\[Beta]=2",pcLabelStyle],Scaled[{0.70,0.92}]],
        Inset[Style["\[Beta]=3/2",pcLabelStyle],Scaled[{0.85,0.92}]],
        Inset[Style["\[Beta]=1",pcLabelStyle],Scaled[{0.92,0.8}]],
        Inset[Style["\[Beta]=0",pcLabelStyle],Scaled[{0.94,0.4}]],
        Inset[Style[StringJoin[
          OptionValue["lPanelLabel"],"\n",
          "\!\(\*
StyleBox[\"n\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"=\",\nFontSlant->\"Italic\"]\)"<>ToString[n/2],"\n",
          "\!\(\*
StyleBox[\"r\",\nFontSlant->\"Italic\"]\)="<>(Round[r,0.01]/.{-0.43->"-3/7",0.->"0"})
        ],pcLabelStyle],Scaled@{0.05,0.84},Scaled@{0,1},Alignment->Left]
      },style2,style
    ]
  ]
]


(* ::Chapter:: *)
(*End*)


End[]


Protect[pDecay,lB,pqDiag];


EndPackage[]
