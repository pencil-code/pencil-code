#!/software/sse/manual/mathematica/12.0.0/bin/wolfram -script
(* ::Package:: *)

Needs["pc`"]
pcInitialize[]

sim="your/run/directory";
data=readTS[sim,"t","urms"]//Transpose;

Export[sim<>"/figure.eps",ListPlot[data]]
