(* ::Package:: *)

(* :Name: pcUtils` *)
(* :Author: Hongzhe Zhou, Stockholm, 2022*)
(* :Version: 0.1 *)
(* :Summary:
    Here includes utils that are used by other parts of the package.
*)


BeginPackage["pcUtils`","pcRead1D`"]


(* ::Chapter:: *)
(*Usage messages*)


(*messages*)
pcAround::usage="pcAround[l] gives Around[l] if l={l1,l2,...},
or {Around[x],Around[y],...} if l={{x1,y1,...},{x2,y2,...},...}."

pcDivide::usage="pcDivide[l,n:3] partitions l into n pieces, with the
length of each no larger than Length[l]/n. The rest elements in l are
thrown away."

pcDifferences::usage="pcDifferences[l] returns
{(l2-l1)/2, (l2-l1)/2+(l4-l3)/2, ..., (l[n]-l[n-1])/2}."

pcFit::usage="pcFit[data,sp,fact:{1,1,1}] fits the data with
some model specified by sp, prints out the result, and returns the fitted model.
Inputs:
  data: A List of the form {{x1,y1},{x2,y2},...,{xn,yn}.
  sp: Either a String that matches the implemented models, or simply an expression
      using \"x\" as the variable and \"a\"[i] as parameters.
      Example: \"a\"[1]+Log[\"a\"[2],\"x\"]^(\"a\"[3]).
      For the moment allows for up to 3 paramters.
  fact: Optional. A List of length 3. The fitted curve is rescaled: The x coordinates of the
        first and last points are rescaled by a factor of fact[[1]] and fact[[2]],
        respectively, and the whole curve is rescaled by fact[[3]].
Outputs:
  Prints out the fitted model (PowerLaw and Exp models also return a linear curve),
  and returns an Association object, with the following keys:
    \"FittedModel\" -> the resulting FittedModel object
    \"FittedCurve\" -> the fitted and rescaled data
    \"Parameters\" -> the fitted parameters
    \"Errors\"  -> errors of the fitted parameters
    \"PwE\" -> the fitted parameters with errors around them."

pcNumberToTex::usage="pcNumberToTex[x,n:1] first converts a number x into scientific form
with n-digit precision, and then converts it into a Tex-form string.
Example:
  pcNumberToTex[0.12345,2]=\"$1.2\\times 10^{-1}$\"."

pcWriteTexTable::usage="pcWriteTexTable[file,head,content,Options] writes into file a Tex-form table.
Input:
  file: If already exists, will overwrite it.
  head: A List of Strings. Must be in Tex-form.
  content: A nested List. Must be in Tex-form.
Options:
  \"lAlignment\": By default ->\"c\", center-aligned. Can also be \"l\", left-aligned.
                  Otherwise, can be a String like \"lccc\", whose length must match that
                  of the header.
Output:
  The output file reads something like this:
\\begin{tabular}{cc}
\\hline
% head1 & head2\\
\\hline
% tb11 & tb12 \\
% tb21 & tb22 \\
\\hline
\\end{tabular}"

pcCombine::usage="pcCombine[sim,simCont,readF,args] combines data produced by two runs,
sim and simCont, where simCont is supposed to be a run restarted from some intermediate
time of sim. Data from sim after this restarting time will be removed.
Inputs:
  sim, simCont: String. Directories of the two runs.
  readF: The function needed to read the data file(s).
  args: Argument(s) needed by readF.
Outputs:
  Same format as that of readF.
Example:
  pcCombine[sim,simCont,read1D2Scale,\"power_kin.dat\",5] to combine the spectrum data."

pcMergePower::usage="pcMergePower[sim,simCont] merges the data/power*.dat files of sim
and simCont, where simCont is supposed to be a run restarted from some intermediate
time of sim. Data from sim after this restarting time will be removed. The old power*.dat
files will be renamed to power*.dat.bk. If a .bk file exists, it will not be rewritten and
the corresponding data file will be skipped and not merged."

pcPkgDir::usage="pcPkgDir[] opens the directory of this package."


Begin["`Private`"]


(* ::Chapter:: *)
(*Functions*)


pcAround[l_]:=Around[l]/;Depth[l]==2
pcAround[l_]:=Around/@Transpose[l]/;Depth[l]==3

pcDivide[l_,n_Integer:3]:=Partition[l,UpTo@Floor[Length[l]/n]][[1;;n]]

pcDifferences[l_]:=Module[{dl},
  dl=l//Differences;
  {dl[[1]]/2,Mean/@(Partition[dl,2,1]),dl[[-1]]/2}//Flatten
]

pcFit[data_,sp_,fact_List:{1,1,1}]:=Module[
  {llinear,funcx,funcy,model,a,x,tmp,dx,dy,weights,fit,fittedCurve,p,e},
  llinear=False;
  funcx=funcy={Identity,Identity};
  
  dx=data[[;;,1]]/._?NumericQ->0/.x_Around->x["Uncertainty"];
  dy=data[[;;,2]]/._?NumericQ->0/.x_Around->x["Uncertainty"];
  weights=dx^2+dy^2;
  weights=If[MemberQ[weights,x_/;x==0],Automatic,1/weights];
  tmp=data/.Around[xx_,yy_]:>xx;
  
  model=Switch[sp,
    "PowerLaw",llinear=True;funcx=funcy={Log,Exp},
    "PowerLaw+C",a[1]+a[2]*x^a[3],
    "Linear",llinear=True,
    "Exp",llinear=True;funcy={Log,Exp},
    "Exp+C",a[1]+a[2]*Exp[a[3]*x],
    _,sp/.{"x"->x,"a"->a}
  ];
  tmp=Transpose[{funcx[[1]][tmp[[;;,1]]],funcy[[1]][tmp[[;;,2]]]}];
  
  fit=If[llinear,
    LinearModelFit[tmp,x,x,Weights->weights],
    NonlinearModelFit[tmp,model,a/@Range[3],x,Weights->weights]
  ];
  fittedCurve=Module[{xx,yy},
    xx=fact[[1;;2]]*MinMax[funcx[[2]][tmp[[;;,1]]]];
    xx=Subdivide[Sequence@@xx,32];
    yy=fact[[3]]*funcy[[2]][fit["Function"]/@funcx[[1]][xx]];
    Transpose[{xx,yy}]
  ];
  
  Print["Fit result: ",fit["BestFit"]/.x->"x"];
  p=If[llinear,fit["BestFitParameters"],(a/@Range[3])/.fit["BestFitParameters"]];
  e=If[Length[data]>=3,fit["ParameterErrors"],ConstantArray[0,Length[p]]];
  Association[
    "FittedModel"->fit,
    "FittedCurve"->fittedCurve,
    "Parameters"->p,
    "Error"->e,
    "PwE"->Around@@@Transpose[{p,e}]
  ]
]

pkgDir=$InputFileName//DirectoryName;
pcPkgDir[]:=(Run@StringJoin["open ",pkgDir];)


(* ::Section:: *)
(*Output to Tex*)


pcNumberToTex[x_,n_:1]:=StringJoin[
  "$",
  x//ScientificForm[#,n]&//ToString[#,TeXForm]&//StringReplace[#,".\\"->"\\"]&,
  "$"
]

Options[pcWriteTexTable]={"lAlignment"->"c"};
pcWriteTexTable[file_,head_,content_,OptionsPattern[]]:=Module[{a},
  pcWriteTexTable::badalign="The given alignment `1` is not compatible with the header.";
  a=Switch[OptionValue["lAlignment"],
    "c",StringJoin[ConstantArray["c",Length[head]]],
    "l",StringJoin[ConstantArray["l",Length[head]]],
    _,OptionValue["lAlignment"]];
  If[StringLength[a]!=Length[head],
    Message[pcWriteTexTable::badalign,a];Return[$Failed]];
  
  Print["Overwriting ",file];
  If[FileExistsQ[file],DeleteFile[file]];
  CreateFile[file];
  
  (*start table*)
  WriteString[file,"\\begin{tabular}{"];
  WriteString[file,a];
  WriteString[file,"}\n\\hline\n"];
  
  (*headings line*)
  WriteString[file,StringRiffle[head," & "]];
  WriteString[file,"\\\\\n\\hline\n"];
  
  (*content*)
  Do[
    WriteString[file,StringRiffle[ToString/@str," & "]];
    WriteString[file,"\\\\\n"],
    {str,content}];
  WriteString[file,"\\hline\n"];
  
  (*finish*)
  WriteString[file,"\\end{tabular}"];
  FilePrint[file]
]


(* ::Section:: *)
(*Combine data*)


pcCombine[sim_String,simCont_String:"None",readF_,args__]:=
Module[{sim2,data1,data2,t0,pos},
  pcCombine::noNewDir="The second directory is not given, and newDirCont.in is not found.";
  
  sim2=simCont;
  If[sim2=="None",
    If[FileExistsQ[sim<>"/newDirCont.in"],
      sim2=Import[sim<>"/newDirCont.in"],
      Message[pcCombine:noNewDir];Return[$Failed]
    ]
  ];
  
  data1=readF[sim,args];
  data2=readF[sim2,args];
  
  (* remove the t>=t0 part in data1 *)
  t0=data2[[1,1]];
  pos=Position[data1[[1]],tt_/;tt<t0];
  data1=Extract[#,pos]&/@data1;
  
  Join@@@Transpose[{data1,data2}]
]

pcMergePower[sim_String,simCont_String]:=
Module[{files,old,bk,combined},
  pcMergePower::oldExists="A back-up file .old exists for `1`. Skipping.";
  
  files=FileNameTake/@FileNames["power*",sim<>"/data/"];
  
  Do[
    old=FileNameJoin[{sim,"data",ff}];
    bk=FileNameJoin[{sim,"data",ff<>".old"}];
    
    If[FileExistsQ[bk],Message[pcMergePower::oldExists,ff];Continue[]];
    combined=pcCombine[sim,simCont,read1D,ff];
    RenameFile[old,bk];
    Export[old,Riffle@@combined],
    
    {ff,files}
  ];
  
  Print["Finished combining: ",files]
]


(* ::Chapter:: *)
(*End*)


End[]


Protect[
  pcAround,pcDivide,pcDifferences,pcFit,
  pkgDir,
  pcNumberToTex,pcWriteTexTable,
  pcCombine, pcMergePower
]


EndPackage[]
