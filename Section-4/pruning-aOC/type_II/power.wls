#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Title:: *)
(*Compute the percolation fraction in ''AND'' out-component *)


(* ::Text:: *)
(*------ How to run it------ -*)
(*  WolframKernel  - script cavity.wls*)
(*This script generates the heatmap with the largest solution of the cavity theory*)
(*---Code-- -*)


(* ::Chapter:: *)
(*in-degree of TF is fat-tail distribution*)


P[k_,\[Lambda]_]:=(1/k^\[Lambda]-1/(k+1)^\[Lambda])


(* ::Chapter:: *)
(*in-degree of genes is 1+ Poisson*)


SetDirectory[NotebookDirectory[]]
d=.
c=.
\[Alpha]=.
p=.
G[x_,d_,p_]=p \!\(
\*SubsuperscriptBox[\(\[Sum]\), \(k = 1\), \(\[Infinity]\)]\(PDF[PoissonDistribution[d], k - 1] \((1 - 
\*SuperscriptBox[\((1 - x)\), \(k\)])\)\)\);
TF[g_,\[Lambda]_]= \!\(
\*SubsuperscriptBox[\(\[Sum]\), \(k = 1\), \(\[Infinity]\)]\(P[k, \[Lambda]] \((
\*SuperscriptBox[\(g\), \(k\)])\)\)\);


d=0.4;
p=0.95;
\[Lambda]min= 1.8;
\[Lambda]max = 5;
\[Lambda]step= 0.1;
list={d,Table[{\[Lambda],x/.NSolve[G[TF[x,\[Lambda]],d,p]-x==0&&x>=0&&x<=1,x,Reals][[-1]]},{\[Lambda],\[Lambda]min,\[Lambda]max,\[Lambda]step}]}
Export["./genes_on_cavity.txt",list,"Table"]


(* ::Text:: *)
(*compute the fraction of TF in aoC*)


list={d,Table[{\[Lambda],x/.NSolve[TF[G[x,d,p],\[Lambda]]-x==0&&x>=0&&x<=1,x,Reals][[-1]]},{\[Lambda],\[Lambda]min,\[Lambda]max,\[Lambda]step}]}
Export["./TFs_on_cavity.txt",list,"Table"]









