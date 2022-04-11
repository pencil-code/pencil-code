(* ::Package:: *)

(* :Name: pcGetParam` *)
(* :Author: Hongzhe Zhou, Stockholm, 2021*)
(* :Version: 0.1 *)
(* :Summary:
    Defines getParam[] which produces various useful parameters.
*)


BeginPackage["pcGetParam`","pcReadBasic`","pcRead1D`"]


(* ::Chapter:: *)
(*Usage messages*)


pcReload::usage="Reload all pc`variables that use Set[].";

getParam::usage="getParam[sim,var] gives various parameters of sim.
Input:
  sim: String. Directory of the simulation folder
  var: String. Specifies the needed variable
Example:
  getParam[dir,\"ReN\"] returns the Reynolds number.";

LuNspec::usage="LuNspec[sim] computes time-independent Lundquist numbers from power_mag.dat.
The definition B*l^(n-1)/eta, where:
B is the rms magnetic energy exlucding the k=0 mode;
The length scale l is the energy-weighted value of 1/k;
n is the power in the diffusion operator del^n (normal diffusion n=2).
Input:
  sim: String. Directory of the simulation folder
Output:
  A time series of the Lundquist number."

paramMatrix::usage="paramMatrix[sim,{file1,var1},{file2,var2},...] gives in the 
MatrixForm values of {var1,var2,...} that can be retrived by readParamNml."

infoMatrix::usage="infoMatrix[sim,var1,var2,...] gives in the MatrixForm values of
{var1,var2,...} that can be retrieved by getParam[sim,\"var\"].";


Begin["`Private`"]


(* ::Chapter:: *)
(*Functions*)


(* ::Section:: *)
(*Private variables (Set[] is only allowed here!)*)


(*some frequently read variables; better to dynamically define*)
(*when data file is updated, pcReload[] is needed to clear existing definitions*)
pcReload[]:=Module[{},
  Clear[urms,urmskf];
  urms[sim_]:=urms[sim]=With[{u=readTS[sim,"urms"]},u[[-Round[Length[u]/3];;-1]]//Mean];
  (*urms of modes with k\[GreaterEqual]kf*)
  urmskf[sim_]:=urmskf[sim]=With[{file=sim<>"/data/power_kin.dat"},
    If[FileExistsQ[file],
      Sqrt[2*Total[(read1D2Scale[sim,"power_kin.dat"])//Last//Last]],
      "NaN"
    ]];
  
  Clear[hrms];
  hrms[sim_]:=hrms[sim]=With[{h=readTS[sim,"oum"]},h[[-Round[Length[h]/3];;-1]]//RootMeanSquare];
  
  Clear[kf];
  kf[sim_]:=kf[sim]=
    If[readParamNml[sim,"start.in","lforcing"] && readParamNml[sim,"run.in","FORCE"]!=0.,
      If[Length[#[[1]]]==1,#[[2,1]],#[[1,2]]]&@Import[FileNameJoin[{sim,"k.dat"}]],
      If[MemberQ[readParamNml[sim,"run.in","IFORCING_CONT"],x_/;StringStartsQ[x,"'GP_TC13"]],
        readParamNml[sim,"run.in","KMIN_GP"],
        "No forcing"
      ]
    ];
  
  Clear[ted,tedkf];
  ted[sim_]:=ted[sim]=1/(kf[sim]*urms[sim]);	
  tedkf[sim_]:=tedkf[sim]=1/(kf[sim]*urmskf[sim]);
  
  Clear[nu];
  nu[sim_]:=nu[sim]=readParamNml[sim,"run.in","NU"];
  
  Clear[eta,etaTFM];
  eta[sim_]:=eta[sim]=
    If[readParamNml[sim,"start.in","lmagnetic"],
      readParamNml[sim,"run.in","ETA"],
      "No magnetic"
    ];
  etaTFM[sim_]:=readParamNml[sim,"run.in","ETATEST"];
  
  Clear[kappaTFM];
  kappaTFM[sim_]:=readParamNml[sim,"run.in","KAPPATEST"];
  
  Clear[PrM];
  PrM[sim_]:=PrM[sim]=nu[sim]/eta[sim];
  
  Clear[omega];
  omega[sim_]:=omega[sim]=Module[{os,or},
    os=readParamNml[sim,"start.in","OMEGA"];
    or=readParamNml[sim,"run.in","OMEGA"];
    Which[
      os==or, Return[os],
      os==0 && or!=0, Return[or],
      os!=0 && or==0, Return[os],
      os!=0 && or!=0 && os!=or, 
        Print["Warning: different OMEGA in start.in and run.in; returning the latter."]];Return[or]
    ]
]
pcReload[]


(* ::Section:: *)
(*getParam*)


(*directly readables*)
getParam[sim_,"folder"]:=StringJoin["...",StringTake[sim,-7]]
getParam[sim_,"id"]:=
  Module[{f},
    If[And[FileExistsQ[FileNameJoin[{sim,"sim_comment"}]],(f=Import[FileNameJoin[{sim,"sim_comment"}]])=!=""],
      First@If[Length[f]==2,Flatten/@f,StringSplit[f,"\n"]],
      "NaN"
    ]
  ]				
getParam[sim_,"comment"]:=
  Module[{f},
    If[And[FileExistsQ[FileNameJoin[{sim,"sim_comment"}]],(f=Import[FileNameJoin[{sim,"sim_comment"}]])=!=""],
      Rest@If[Length[f]==2,Flatten/@f,StringSplit[f,"\n"]],
      "No comment."
    ]
  ]

getParam[sim_,"Nx"]:=readDim[sim]["nx"]
getParam[sim_,"Ny"]:=readDim[sim]["ny"]
getParam[sim_,"Nz"]:=readDim[sim]["nz"]
getParam[sim_,"Nxyz"]:=readDim[sim]/@{"nx","ny","nz"}
getParam[sim_,"LXYZ"]:=readParamNml[sim,"start.in","LXYZ"]

getParam[sim_,"urms"]:=urms[sim]
getParam[sim_,"urmskf"]:=urmskf[sim]
getParam[sim_,"hrms"]:=hrms[sim]
getParam[sim_,"kf"]:=kf[sim]
getParam[sim_,"nu"]:=nu[sim]
getParam[sim_,"eta"]:=eta[sim]

(*dimensionless numbers*)
getParam[sim_,"ReN"]:=urms[sim]/kf[sim]/nu[sim]
getParam[sim_,"ReN",k2_]:=urms[sim]/k2/nu[sim] (*supply kf by hand*)
getParam[sim_,"ReNkf"]:=urmskf[sim]/kf[sim]/nu[sim]

getParam[sim_,"ReM"]:=urms[sim]/kf[sim]/eta[sim]
getParam[sim_,"ReM",k2_]:=urms[sim]/k2/eta[sim] (*supply kf by hand*)
getParam[sim_,"PrM"]:=PrM[sim]
getParam[sim_,"PrMTFM"]:=nu[sim]/etaTFM[sim]

getParam[sim_,"Ro"]:=If[omega[sim]==0,"No rotation",kf[sim]*urms[sim]/2/omega[sim]]
getParam[sim_,"Ro",k2_]:=If[omega[sim]==0,"No rotation",k2*urms[sim]/2/omega[sim]]
getParam[sim_,"Rokf"]:=If[omega[sim]==0,"No rotation",kf[sim]*urmskf[sim]/2/omega[sim]]

getParam[sim_,"Co"]:=If[omega[sim]==0,0,1/getParam[sim,"Ro"]]
getParam[sim_,"Co",k2_]:=If[omega[sim]=0,0,1/getParam[sim,"Ro",k2]]
getParam[sim_,"Cokf"]:=If[omega[sim]==0,0,1/getParam[sim,"Rokf"]]

getParam[sim_,"Sh"]:=
  If[readParamNml[sim,"start.in","lshear"],
    readParamNml[sim,"run.in","SSHEAR"]/urms[sim]/kf[sim],
    0
  ]
getParam[sim_,"Shkf"]:=
  If[readParamNml[sim,"start.in","lshear"],
    readParamNml[sim,"run.in","SSHEAR"]/urmskf[sim]/kf[sim],
    0
  ]
getParam[sim_,"-Sh"]:=-getParam[sim,"Sh"]
getParam[sim_,"-Shkf"]:=-getParam[sim,"Shkf"]
getParam[sim_,"Sh",k2_]:=
  If[readParamNml[sim,"start.in","lshear"],
    readParamNml[sim,"run.in","SSHEAR"]/urms[sim]/k2,
    0
  ]
getParam[sim_,"-Sh",k2_]:=-getParam[sim,"Sh",k2]
getParam[sim_,"qshear"]:=
  Module[{s=readParamNml[sim,"run.in","SSHEAR"],o=readParamNml[sim,"start.in","OMEGA"]},
    If[s===$Failed,s=0];
    If[s==0&&o==0,"/"//Return];
    If[s!=0&&o==0,"No rotatoin"//Return];
    Rationalize[-s/o]/.Rational[x_,y_]:>ToString[x]<>"/"<>ToString[y]
]

getParam[sim_,"ScTFM"]:=nu[sim]/kappaTFM[sim]

(*secondaries*)
(*eddy turnover time*)
getParam[sim_,"ted"]:=ted[sim]
getParam[sim_,"ted",k2_]:=2\[Pi]/urms[sim]/k2
getParam[sim_,"tedkf"]:=tedkf[sim]
(*diffusive timescales*)
getParam[sim_,"teta1",k1_:1]:=1/(eta[sim]*k1^2)
(**)
getParam[sim_,"knu"]:=kf[sim]*getParam[sim,"ReN"]^(3/4)
getParam[sim_,"knu",k2_]:=k2*getParam[sim,"ReN",k2]^(3/4) (*supply kf by hand*)
getParam[sim_,"knumax"]:=(Max[readTS[sim,"urms"]]/nu[sim])^(3/4)*kf[sim]^(1/4)
getParam[sim_,"knumax",k2_]:=(Max[readTS[sim,"urms"]]/nu[sim])^(3/4)*k2^(1/4) (*supply kf by hand*)
getParam[sim_,"keta"]:=
  If[eta[sim]=="No magnetic","No magnetic",
    getParam[sim,"knu"]/If[PrM[sim]>=1,PrM[sim]^(-1/2),PrM[sim]^(-3/4)]
  ]
getParam[sim_,"keta",k2_]:=
  If[eta[sim]=="No magnetic","No magnetic",
    getParam[sim,"knu",k2]/If[PrM[sim]>=1,PrM[sim]^(-1/2),PrM[sim]^(-3/4)]
  ]
getParam[sim_,"ketamax"]:=
  If[eta[sim]=="No magnetic","No magnetic",
    getParam[sim,"knumax"]/If[PrM[sim]>=1,PrM[sim]^(-1/2),PrM[sim]^(-3/4)]
  ]
getParam[sim_,"ketamax",k2_]:=
  If[eta[sim]=="No magnetic","No magnetic",
    getParam[sim,"knumax",k2]/If[PrM[sim]>=1,PrM[sim]^(-1/2),PrM[sim]^(-3/4)]
  ]
getParam[sim_,"kRo"]:=If[omega[sim]==0,"No rotation",kf[sim]*getParam[sim,"Ro"]^(-3/2)]
getParam[sim_,"kRo",k2_]:=If[omega[sim]==0,"No rotation",k2*getParam[sim,"Ro"]^(-3/2)]


(* ::Section:: *)
(*Dimensionless parameters from spectra files*)


LuNspec[sim_]:=Module[{t,spec,Eb,k,l,n,eta},
  (*error messages*)
  LuNspec::nofile="power_mag.dat not found from `1`.";
  LuNspec::nores="Unfamiliar iresistivity for `1`.";
  
  (**)
  If[!FileExistsQ[sim<>"/data/power_mag.dat"],
    Message[LuNspec::nofile,sim];Return[$Failed]
  ];
  {t,spec}=read1D[sim,"power_mag.dat"];
  If[t[[1]]==0.,{t,spec}=Rest/@{t,spec}];
  spec=Rest/@spec;
  Eb=2Total/@spec;
  l=(Total[1/Range[spec//First//Length]*#]&/@spec)/Eb;
  
  {eta,n}=Switch[readParamNml[sim,"run.in","IRESISTIVITY"]//First,
    "'eta-const'",  {readParamNml[sim,"run.in","ETA"],       2},
    "'eta-tdep'",   {readParamNml[sim,"run.in","ETA"],       2},
    "'hyper2'",     {readParamNml[sim,"run.in","ETA_HYPER2"],4},
    "'hyper3'",     {readParamNml[sim,"run.in","ETA_HYPER3"],6},
    "'hyper3-tdep'",{readParamNml[sim,"run.in","ETA"],       6},
    _,               Message[LuNspec::nores,sim];Return[$Failed]
  ];
  
  {t,Sqrt[Eb]*l^(n-1)/eta}//Transpose
]


(* ::Section:: *)
(*Handy functions*)


paramMatrix[sim_,vars__]:=Module[{files,prms,data},
  {files,prms}=Transpose@{vars};
  data=MapThread[readParamNml[sim,#1,#2]&,{files,prms}];
  MatrixForm[Transpose@{prms,data}]
]

infoMatrix[sim_,vars__]:=MatrixForm[Transpose[{{vars},(getParam@@Flatten[{sim,#}])&/@{vars}}]]


(* ::Chapter:: *)
(*End*)


End[]


Protect[
  pcReload,getParam,
  LuNspec,
  paramMatrix,infoMatrix
]


EndPackage[]
