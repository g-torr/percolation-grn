#!/usr/local/bin/WolframKernel
(* ::Package:: *)

(* ::Text:: *)
(*#------ How to run it------ -*)
(* WolframKernel  - script cavity.wls*)


(* ::Text:: *)
(*This script generates the heatmap with the largest solution of the cavity theory*)


(* ::Text:: *)
(*---Code-- -*)


SetDirectory[NotebookDirectory[]]
Print["Data will be saved in: ",Directory[], "/on_cavity.txt"]
cmin=0.01;
cmax=6;
cstep=0.05;
dmin=0.01;
dmax=6;
dstep=0.05;
p=0.95;
G[x_,d_,p_]= p \!\(
\*SubsuperscriptBox[\(\[Sum]\), \(k = 1\), \(\[Infinity]\)]\(PDF[PoissonDistribution[d], k - 1] \((1 - 
\*SuperscriptBox[\((1 - x)\), \(k\)])\)\)\);
TF[x_,c_]=\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(k = 1\), \(\[Infinity]\)]\(PDF[PoissonDistribution[c], k - 1] \((
\*SuperscriptBox[\(x\), \(k\)])\)\)\);




list = Table[{d, 
  	 Table[{c, 
    		g /. NSolve[
      				G[TF[g, c], d,  p] - g == 0 && g >= 0 && 
       g <= 1, g, Reals][[-1]]}, {c, cmin, cmax, cstep}]}, {d, dmin, dmax, 
  dstep}]
Export["./genes_on_cavity.txt", list, "Table"]
Print["Doing TFs now"]
G[x_,d_,p_]= p \!\(
\*SubsuperscriptBox[\(\[Sum]\), \(k = 1\), \(\[Infinity]\)]\(PDF[PoissonDistribution[d], k - 1] \((1 - 
\*SuperscriptBox[\((1 - x)\), \(k\)])\)\)\);
TF[x_,c_]=\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(k = 1\), \(\[Infinity]\)]\(PDF[PoissonDistribution[c], k - 1] \((
\*SuperscriptBox[\(x\), \(k\)])\)\)\);
list = Table[{d, 
  	 Table[{c, 
    		t /. NSolve[
      				TF[G[t, d,p], c] - t == 0 && t >= 0 && 
       t <= 1, t, Reals][[-1]]}, {c, cmin, cmax, cstep}]}, {d, dmin, dmax, 
  dstep}]
Export["./TFs_on_cavity.txt", list, "Table"]


Print["Data will be saved in: ",Directory[], "/genes_off_cavity.txt"]
d =.
c =.
p =.



(* ::Text:: *)
(*Genes _initially Off*)


G2[x_, d_,  p_] = (1 - p) \!\(
\*SubsuperscriptBox[\(\[Sum]\), \(k = 1\), \(\[Infinity]\)]\(PDF[PoissonDistribution[d], k - 1] \((1 - 
\*SuperscriptBox[\((1 - x)\), \(k\)])\)\)\) + p
TF2[x_, c_] = \!\(
\*SubsuperscriptBox[\(\[Sum]\), \(k = 1\), \(\[Infinity]\)]\(PDF[PoissonDistribution[c], k - 1] \((
\*SuperscriptBox[\(x\), \(k\)])\)\)\);
p=0.95;
list = Table[{d, 
		Table[{c, 
     			g/.NSolve[G2[TF2[g, c], d,  1-p] - g == 0 && g >= 0 &&g <= 1, g, Reals][[1]]},  {c, cmin, cmax,dstep}]},
		 {d, dmin, dmax,cstep}]

Export["./genes_off_cavity.txt", list, "Table"]



list = Table[{d, 
		Table[{c, 
     			g/.NSolve[TF2[G2[g, d, 1-p], c] - g == 0 && g >= 0 &&g <= 1, g, Reals][[1]]},  {c, cmin, cmax,dstep}]},
		 {d, dmin, dmax,dstep}]

Export["./TFs_off_cavity.txt", list, "Table"]


