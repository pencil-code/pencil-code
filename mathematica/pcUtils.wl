(* ::Package:: *)

(* :Name: pcUtils` *)
(* :Author: Hongzhe Zhou, Stockholm, 2022*)
(* :Version: 0.1 *)
(* :Summary:
    Here includes utils that are used by other parts of the package.
*)


BeginPackage["pcUtils`"]


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

pcFit::usage="pcFit[data,sp,fact:{1,1,1},\"lEcho\"->False] fits the data with 
some model specified by sp, prints out the result, and returns a the fitted curve.
One can also use Reap[pcFit[...]] to get the fitting parameters.
Inputs:
  data: A List of the form {{x1,y1},{x2,y2},...,{xn,yn}.
  sp: Either a String that matches the implemented models, or simply an expression
      using \"x\" as the variable and \"a\"[i] as parameters.
      Example: \"a\"[1]+Log[\"a\"[2],\"x\"]^(\"a\"[3]).
      For the moment allows for up to 10 paramters.
  fact: Optional. A List of length 3. The fitted curve is rescaled: The x coordinates of the
        first and last points are rescaled by a factor of fact[[1]] and fact[[2]],
        respectively, and the whole curve is rescaled by fact[[3]].
Options:
  \"lEcho\": If True, then returns the original data rather than the fitted curve.
Outputs:
  Prints out the fitted parameters, and returns a List."

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

Options[pcFit]={"lEcho"->False};
pcFit[data_,sp_,fact_List:{1,1,1},OptionsPattern[]]:=Module[{model,llog,a,x,sol,minmax},
  llog=False;
  model=Switch[sp,
    "PowerLaw",llog=True;a[1]+a[2]*x,
    "PowerLaw+C",a[1]+a[2]*x^a[3],
    "Linear",a[1]+a[2]*x,
    "Exp",a[1]*Exp[a[2]*x],
    "Exp+C",a[1]+a[2]*Exp[a[3]*x],
    _,sp/.{"x"->x,"a"->a}
  ];
  Sow[model/.{a[i_]:>("a["<>ToString@i<>"]"),x->"x"}];
  If[llog,
    sol=FindFit[Log@data,model,a/@Range[10],x];
    model=Exp[model/.x->Log[x]/.sol],
    (*else*)
    sol=FindFit[Log@data,model,a/@Range[10],x];
    model=model/.sol
  ];
  Sow[sol/.a[i_]:>("a["<>ToString@i<>"]")];
  Print["Fit result: ",model/.x->"x"];
  If[OptionValue["lEcho"],data//Return];
  minmax=MinMax[data[[;;,1]]]*fact[[1;;2]];
  Table[{x,fact[[3]]*model},{x,Subdivide[Sequence@@minmax,32]}]
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


(* ::Chapter:: *)
(*End*)


End[]


Protect[
  pcAround,pcDivide,pcDifferences,pcFit,
  pkgDir,
  pcNumberToTex,pcWriteTexTable
]


EndPackage[]
