(* ::Package:: *)

(********************************************************************)
(*                                                                  *)
(* :Title:       GridMethod2D_PA                                    *)
(* :Authors:     Chaney Lin                                         *)
(* :Version:     Unofficial                                         *)
(* :Date:        May 2015                                           *)
(*                                                                  *)
(********************************************************************)

BeginPackage["Tilings`GridMethod2DPA`"]

(********************************************************************)
(*              Usage messages for exported functions               *)
(********************************************************************)

StarVectors::usage=
"
"

(********************************************************************)
(*                                                                  *)
(*                   Unprotect exported functions                   *)
(*                                                                  *)
(********************************************************************)

Unprotect[GetEffectiveN,
	GetMinDist,
	StarVectors,
	PlotGridVectors,
	PlotGrid,
	GetIndexList,
	GetIntersections,
	PlotIntersections,
	GetIntersectionTiles,
	DualizeGrid,
	GetSingleTile,
	PlotDualTiling,
	PlotSingular,
	PlotColors,
	PlotColors2,
	GetPoints,
	GetPoints2,
	GetLines,
	GetTiles,
	GetMainVertex,
	GetLines2,
	GetTiles2,
	GetMainVertex2,
	PlotIntersectionTiles,
	GetUnitCell,
	PlotUnitCell,
	add,
	add2,
	TrimUnitCell,
	GetAmmannUnitCell,
	FindNeighbors,
	PlotStarVectors,
	StarPhases2,
	StarVectors2,
	StarVectors3,
	DualizeGrid2,
	PlotDualTiling2,
	PlotGrid2,
	PlotIntersections2,
	GetUnitCell2,
	TrimUnitCell2,
	TrimUnitCell3,
	PlotUnitCell2,
	PlotUnitCell3,
	GetUnitCell3,
	GetKLattice,
	DualizeUnitCell,
	PlotDualUnitCell,
	PlotBZ,
	GetInfoGraphics,
	ModStarVectors]


(********************************************************************)
(*                                                                  *)
(*                     Start of Private context                     *)
(*                                                                  *)
(********************************************************************)

Begin["`Private`"]

Needs["ComputationalGeometry`"];

(*********************************)
(* Clearing previous definitions *)
(*********************************)

Clear[GetEffectiveN,
	GetMinDist,
	StarVectors,
	PlotGridVectors,
	PlotGrid,
	GetIntersections,
	GetIndexList,
	PlotIntersections,
	GetIntersectionTiles,
	DualizeGrid,
	GetSingleTile,
	PlotDualTiling,
	PlotSingular,
	PlotColors,
	PlotColors2,
	GetPoints,
	GetPoints2,
	GetLines,
	GetTiles,
	GetMainVertex,
	GetLines2,
	GetTiles2,
	GetMainVertex2,
	PlotIntersectionTiles,
	GetUnitCell,
	PlotUnitCell,
	add,
	add2,
	TrimUnitCell,
	GetAmmannUnitCell,
	FindNeighbors,
	PlotStarVectors,
	StarPhases2,
	StarVectors2,
	DualizeGrid2,
	PlotDualTiling2,
	PlotGrid2,
	PlotIntersections2,
	GetUnitCell2,
	TrimUnitCell2,
	TrimUnitCell3,
	PlotUnitCell2,
	PlotUnitCell3,
	GetUnitCell3,
	GetKLattice,
	DualizeUnitCell,
	PlotDualUnitCell,
	PlotBZ,
	GetInfoGraphics,
	StarVectors3,
	ModStarVectors]



$MinPrecision = $MachinePrecision;

add[p_,r_]:=Apply[Plus,p*r];

GetIndexList[theta_,angles_,acc_:.05,out_:False]:=
	Module[{r={},indexset={},curangle=theta,A,invA,alpha,alphaPA,rapprox,rmag,
			inda1,inda2,a1scale,a2scale},

			If[theta+Total[angles]>360,Return[{}];];
			
			AppendTo[r,{1,0}];
			AppendTo[r,{Cos[theta/360*2*Pi],Sin[theta/360*2*Pi]}];

			Do[curangle=a+curangle;
				AppendTo[r,{Cos[curangle/360*2*Pi],Sin[curangle/360*2*Pi]}];
			,{a,angles}];

			A=Transpose[{r[[1]],r[[2]]}];
			invA=Inverse[A];

			Do[
				alpha=invA.ri;
				alphaPA=Rationalize[#,acc]&/@alpha;
				AppendTo[indexset,alphaPA];
			,{ri,r}
			];
			rapprox=A.#&/@indexset;
			rmag=Norm/@rapprox;

			inda1=DeleteCases[indexset[[All,2]],0];
			inda2=DeleteCases[indexset[[All,1]],0];

			a1scale=LCM[Sequence@@(1/#&/@inda1)];
			a2scale=LCM[Sequence@@(1/#&/@inda2)];

			If[out,
				Print["(x,y):",MatrixForm[Transpose[{indexset[[All,1]],indexset[[All,2]]}]]];
				Print["scale:",{a1scale,a2scale}];
			];
			
			a1scale = a1scale / Sin[theta/360*2*Pi];
			a2scale = a2scale / Sin[theta/360*2*Pi];

			{indexset[[{3,4,5}]],{a1scale,a2scale}}
]


GetIndexListFib[theta_,coeffs_,n_,out:False]:=
	Module[{NN,fiblist,indexset,inda1,inda2,a1scale,a2scale},
		
			NN = Length[coeffs[[1,1]]];
		
			fiblist = Table[Fibonacci[j],{j,Range[n+1,NN+n]}];
		
			indexset = Map[fiblist.#&,coeffs,{2}] / Fibonacci[NN + n + 1];

			inda1=DeleteCases[indexset[[All,2]],0];
			inda2=DeleteCases[indexset[[All,1]],0];

			a1scale=LCM[Sequence@@(1/#&/@inda1)];
			a2scale=LCM[Sequence@@(1/#&/@inda2)];
		
			a1scale = a1scale / Sin[theta/360*2*Pi];
			a2scale = a2scale / Sin[theta/360*2*Pi];

	{indexset,{a1scale,a2scale}}
]


StarVectors3[theta_,coeffs_,n_,out_:False]:=
	Module[{indices,scale,rr1,rr2,r1,r2,r3,r4,r5,a1,a2,a3,a4,a5,r,a,rmag,rorth,rangles={},
			NN,tau,tautable,MapCoeffs,Pairs,x1,x2,y1,y2,dot,det,angle},

			{indices,scale}=GetIndexListFib[theta,coeffs,n,out];

			If[indices=={},Return[{}];];

			rr1={1,0};
			rr2={Cos[theta/360*2*Pi],Sin[theta/360*2*Pi]};
			r1=indices[[1,1]]*rr1+indices[[1,2]]*rr2;
			r2=indices[[2,1]]*rr1+indices[[2,2]]*rr2;
			r3=indices[[3,1]]*rr1+indices[[3,2]]*rr2;
			r4=indices[[4,1]]*rr1+indices[[4,2]]*rr2;
			r5=indices[[5,1]]*rr1+indices[[5,2]]*rr2;
			r={r1,r2,r3,r4,r5};
		
			rmag = Map[Norm,r];
			r = Map[Normalize,r];
			r = Table[r[[i]]/rmag[[i]],{i,Range[5]}];
			rmag = Map[Norm,r];
			rorth=Table[{-x[[2]],x[[1]]},{x,r}];

			NN = Length[coeffs[[1,1]]];
			tau=GoldenRatio;
			tautable=N[Table[tau^j,{j,Range[NN]}]/tau^(NN+1)];
			MapCoeffs=Map[tautable.#&,coeffs,{2}];

			a={a1,a2,a3,a4,a5}=#.{rr1,rr2}&/@MapCoeffs;

			Pairs={{1,2},{2,3},{3,4},{4,5},{5,1}};
			
			Do[
				{x1,y1}=r[[p[[1]]]];
				{x2,y2}=r[[p[[2]]]];
				dot=x1*x2+y1*y2;
				det=x1*y2-y1*x2;
				angle=Mod[ArcTan[dot,det],2Pi]/Degree;
				AppendTo[rangles,angle];
			,{p,Pairs}
			];
			
	{N[{r,rorth,rmag,rangles}],N[a],{indices,N[scale]}}
]


StarVectors2[theta_,angles_,acc_:.05,out_:False]:=
	Module[{indices,scale,r1,r2,r3,r4,r5,r,rmag,rorth,rangles={},
			Pairs,x1,y1,x2,y2,dot,det,angle,a1,a2,a3,a4,a5,a},

			{indices,scale}=GetIndexList[theta,angles,acc,out];

			If[indices=={},Return[{}];];

			r1=a1={1,0};
			r2=a2={Cos[theta/360*2*Pi],Sin[theta/360*2*Pi]};
			r3=indices[[1,1]]*r1+indices[[1,2]]*r2;
			r4=indices[[2,1]]*r1+indices[[2,2]]*r2;
			r5=indices[[3,1]]*r1+indices[[3,2]]*r2;
			r={r1,r2,r3,r4,r5};
		
			rmag = Map[Norm,r];
			r = Map[Normalize,r];
			r = Table[r[[i]]/rmag[[i]],{i,Range[5]}];
			rmag = Map[Norm,r];
			rorth=Table[{-x[[2]],x[[1]]},{x,r}];
			
			a3=RotationMatrix[angles[[1]]/360*2*Pi].a2;
			a4=RotationMatrix[angles[[2]]/360*2*Pi].a3;
			a5=RotationMatrix[angles[[3]]/360*2*Pi].a4;
			a={r1,r2,a3,a4,a5};

			Pairs={{1,2},{2,3},{3,4},{4,5},{5,1}};
			
			Do[
				{x1,y1}=r[[p[[1]]]];
				{x2,y2}=r[[p[[2]]]];
				dot=x1*x2+y1*y2;
				det=x1*y2-y1*x2;
				angle=Mod[ArcTan[dot,det],2Pi]/Degree;
				AppendTo[rangles,angle];
			,{p,Pairs}
			];
			
			{N[{r,rorth,rmag,rangles}],N[a],{indices,N[scale]}}
]

PlotStarVectors[starvectors_,modvectors_]:=
	Module[{r,a,vecgfx,origgfx,modgfx},
		r=starvectors[[1,1]];
		a=starvectors[[2]];
		vecgfx=Graphics[{Red,Arrow[{{0,0},#}]&/@r}];
		origgfx=Graphics[Arrow[{{0,0},#}]&/@a];
		modgfx=Graphics[{Blue,Arrow[{{0,0},#}]&/@modvectors}];
		Show[{vecgfx,origgfx,modgfx},PlotRange->{{-1,1},{-1,1}}]
]


GetKListAndGamma[kmin_,kmax_,gamma_:"random",nj_:5,gammaseed_:1]:=
	Module[{kmn,kmx,gam},
		SeedRandom[gammaseed];
		kmn       = Which[And[Head[kmin]===List,
                             Length[kmin]===nj,
                             Union[Map[IntegerQ,kmin]]==={True}],
                         kmin,
                         IntegerQ[kmin],
                         Table[kmin,{nj}],
                         True,
                         Print["Error: invalid argument kmin"];
                         Return[]];
		kmx       = Which[And[Head[kmax]===List,
                             Length[kmax]===nj,
                             Union[Map[IntegerQ,kmax]]==={True}],
                         kmax,
                         IntegerQ[kmax],
                         Table[kmax,{nj}],
                         True,
                         Print["Error: invalid argument kmax"];
                         Return[]];
		gam       = Which[gamma==="random",
                         Print["note: translation vector gamma ",
                               "is chosen randomly!"];
                         Table[Random[Real],{nj}],
                         And[Head[gamma]===List,
                             Length[gamma]==2,
                             First[gamma]==="random",
                             NumberQ[Last[gamma]]],
                         Apply[Append[#,gamma[[2]]-Apply[Plus,#]]&,
                               {Table[Random[Real,{-0.5,0.5}],{nj-1}]}],
                         And[Head[gamma]===List,
                             Length[gamma]===nj,
                             Union[Map[NumberQ,gamma]]==={True}],
                         gamma,
                         NumberQ[gamma],
                         Table[gamma,{nj}],
                         True,
                         Print["Error: invalid argument gamma"];
                         Return[]];
	{kmn,kmx,gam}
]

PlotGrid[kmin_,kmax_,starvectors_,gamma_:"random",gammaseed_,len_:20]:= 
	Module[{r,rorth,rmag,kmn,kmx,gam,g,nj=5},
		{r,rorth,rmag} = starvectors[[1,{1,2,3}]];
		{kmn,kmx,gam}=GetKListAndGamma[kmin,kmax,gamma,nj,gammaseed];
	
		g = Table[{(kj-gam[[j]]/rmag[[j]])*r[[j]] - len*rorth[[j]],
                  (kj-gam[[j]]/rmag[[j]])*r[[j]] + len*rorth[[j]]},
                 {j,nj},
                 {kj,kmn[[j]],kmx[[j]]}];
		g = Map[Line,Flatten[g,1]];
	g
]


DualizeGrid2[kmin_,kmax_,starvectors_,
            gamma_:"random",gammaseed_:1] :=
	Module[{nj=5,r,rorth,rmag,kmn,kmx,gam,step,cut,gi,ri,rior,rmi,kmni,kmxi,gj,rj,rjor,rmj,kmnj,kmxj,li,lj,intersec,tij,
			intersecpts={},tilesij={},tiles={},tilpoints,tillines},

		{r,rorth,rmag} = starvectors[[1,{1,2,3}]];
		{kmn,kmx,gam}=GetKListAndGamma[kmin,kmax,gamma,nj,gammaseed];
          step[i_]  := step[i] = Table[If[j==i,1,0],{j,nj}];
          cut[p_]   := MapThread[Max,{MapThread[Min,{p,kmx+1}],kmn}];

		Print[gam];

          Do[
             gi   = gam[[i]];
             ri   = r[[i]]; 
             rior = rorth[[i]];
			 rmi = rmag[[i]];
             kmni = kmn[[i]];
             kmxi = kmx[[i]];
             Do[
                gj   = gam[[j]];
                rj   = r[[j]];
				rjor = rorth[[j]];
				rmj = rmag[[j]];
                kmnj = kmn[[j]];
                kmxj = kmx[[j]];
                Do[
					li = rmi(rmi*ki - gi);
					lj = rmj(rmj*kj - gj);
					
					intersec = 1/Dot[rior,rj]*(lj*rior-li*rjor);
				
                   tij = Table[Ceiling[N[(Dot[intersec,r[[ii]]/rmag[[ii]]]+
                                         gam[[ii]])/rmag[[ii]]]],{ii,nj}];
                   tij[[i]] = ki;
                   tij[[j]] = kj;

					AppendTo[intersecpts,intersec];
					AppendTo[tilesij,{i,j}];
                   AppendTo[tiles,Map[cut,{tij,tij+step[i],
                                           tij+step[i]+step[j],
                                           tij+step[j]}]];
                   ,{ki,kmni,kmxi}, 
                   {kj,kmnj,kmxj}],
                {j,i+1,nj}], 
             {i,nj-1}];
          tilpoints = Union[Flatten[tiles,1]];
          tillines  = Flatten[Map[Partition[#,2,1]&,
                                  Map[Append[#,#[[1]]]&,
                                      tiles]],1];
          tillines  = Union[Map[Sort,tillines]];
          {nj,tilpoints,tillines,tiles,tilesij,intersecpts,kmn,kmx,gam}
]


GetEffectiveN[til_,uc_,mv_,maxdist_]:=
	Module[{pairs,uniquepairs,vectorpairs,diagtotaldiffdist,shortpos,
			countduplicates,CountTiles},
			pairs=til[[5]][[uc[[1]]]];
			uniquepairs=DeleteDuplicates[pairs];
			vectorpairs=Map[mv[[#]]&,uniquepairs];
			diagtotaldiffdist=Map[Norm/@{Total[#],Differences[#]}&,vectorpairs,{1}];
			Print[diagtotaldiffdist];
			shortpos=Position[diagtotaldiffdist,_?(#<maxdist&)];
			
			CountTiles[ijlist_,pairs_]:=
				Module[{count=0,ijcount},
					Do[
						ijcount=Length[Flatten[Position[pairs,ij]]];
						Print[ij," ",ijcount];
						count=count+ijcount;
						,{ij,ijlist}
					];
				count
			];
			countduplicates=CountTiles[uniquepairs[[shortpos[[All,1]]]],pairs];

	Length[uc[[1]]]-(countduplicates)
]


GetMinDist[mv_]:=
	Module[{ijlist,vectorpairs,diagtotaldiffdist},
			ijlist=DeleteDuplicates[Permutations[Range[5],{2}],Sort[#1]==Sort[#2]&];
			vectorpairs=Map[mv[[#]]&,ijlist];
			diagtotaldiffdist=Map[Norm/@{Total[#],Differences[#]}&,vectorpairs,{1}];
	Min[diagtotaldiffdist]/2
]



GetSingleTile[til_,v1_,v2_]:=
Module[{nj,rows,tiles,tilesij,tilpoints,tillines,intersecpts},
			nj=til[[1]];
			rows=Join[Position[til[[5]],{v1,v2}],Position[til[[5]],{v2,v1}]][[All,1]];
			tiles = Part[til[[4]],rows];
			tilesij=Part[til[[5]],rows];
			intersecpts=Part[til[[6]],rows];
		  tilpoints = Union[Flatten[tiles,1]];
          tillines  = Flatten[Map[Partition[#,2,1]&,
                                  Map[Append[#,#[[1]]]&,
                                      tiles]],1];
          tillines  = Union[Map[Sort,tillines]];

		{nj,tilpoints,tillines,tiles,tilesij,intersecpts}
]

PlotIntersections[til_,kmin_,kmax_,starvectors_,gamma_,
			linewidth_:1/500,
			ptsize_:1/100,
			len_:20,
			plotrange_:6] :=
Module[{g,g12,g13,g14,g15,g23,g24,g25,g34,g35,g45,
		points1,points2,points3,points4,points,
		gam,r,rorth,rmag,indices,scale,a1,a2,unitcellvert,center,unitcell},

		g = PlotGrid[kmin,kmax,starvectors,gamma];

		g12=GetSingleTile[til,1,2];
		g13=GetSingleTile[til,1,3];
		g14=GetSingleTile[til,1,4];
		g15=GetSingleTile[til,1,5];
		g23=GetSingleTile[til,2,3];
		g24=GetSingleTile[til,2,4];
		g25=GetSingleTile[til,2,5];
		g34=GetSingleTile[til,3,4];
		g35=GetSingleTile[til,3,5];
		g45=GetSingleTile[til,4,5];

		points1 = Table[Join[{PointSize[ptsize],Green,Opacity[0.8]},Map[Point,x[[6]]]],{x,{g12,g13,g14,g15}}];
		points2 = Table[Join[{PointSize[ptsize],Blue,Opacity[0.8]},Map[Point,x[[6]]]],{x,{g23,g24,g25}}];
		points3 = Table[Join[{PointSize[ptsize],Red,Opacity[0.8]},Map[Point,x[[6]]]],{x,{g34,g35}}];
		points4 = Table[Join[{PointSize[ptsize],Black,Opacity[0.8]},Map[Point,x[[6]]]],{x,{g45}}];

		points={points1,points2,points3,points4};

		gam=til[[9]];
		{r,rorth,rmag}=starvectors[[1,{1,2,3}]];
		{indices,scale}=starvectors[[3]];
		a1=rorth[[1]]*scale[[1]];
		a2=rorth[[2]]*scale[[2]];
		unitcellvert=N[Table[x-a1/2-a2/2+gam[[1]]*Normalize[r[[1]]]+gam[[2]]*Normalize[r[[2]]],{x,{{0,0},a1,a1+a2,a2}}]];
		center=Total[unitcellvert]/4;
		unitcell=Polygon[unitcellvert];

		Show[{Graphics[Join[g,points]],Graphics[{Opacity[0.1],unitcell}]},
			AspectRatio -> Automatic,
			PlotRange->plotrange{{-1,1},{-1,1}}+center]

]


ModStarVectors[angles_]:=
Module[{r={{1,0}},theta},
		theta=0;
		Do[
			theta=theta+angle;
			AppendTo[r,{Cos[theta/360*2*Pi],Sin[theta/360*2*Pi]}];
		,{angle,angles}
		];
	N[r]
]


GetUnitCell2[theta_,til_,starvectors_,modvectors_]:=
Module[{rr1,rr2,rr1o,rr2o,r,rorth,rmag,a,indices,scale,nj,gam,a1,a2,unitcellvert,unitcellvertindices,mapvert,basisvectors,n,unitcellhull,tiles,tilpoints},
		
		{r,rorth,rmag}=starvectors[[1,{1,2,3}]];
		a=starvectors[[2]];
		{indices,scale}=starvectors[[3]];

		
		rr1={1,0};
		rr2={Cos[theta/360*2*Pi],Sin[theta/360*2*Pi]};
		{rr1o,rr2o}=Table[{-x[[2]],x[[1]]},{x,{rr1,rr2}}];

		nj = til[[1]];
		gam=til[[9]];

		a1=rr1o*scale[[1]];
		a2=rr2o*scale[[2]];
		unitcellvert=N[Table[x-a1/2-a2/2+gam[[1]]*Normalize[r[[1]]]+gam[[2]]*Normalize[r[[2]]],{x,{{0,0},a1,a1+a2,a2}}]];
		unitcellvertindices = Table[Table[Ceiling[N[(Dot[x,r[[ii]]/rmag[[ii]]]+
                                         gam[[ii]])/rmag[[ii]]]],{ii,nj}],{x,unitcellvert}];

		mapvert = Map[add[#,modvectors]&,unitcellvertindices];
		basisvectors={mapvert[[3]]-mapvert[[4]], mapvert[[1]]-mapvert[[4]]};		
		
		unitcellhull = ConvexHull[unitcellvert];
		indices={};
		n = Length[til[[6]]];
		Do[If[unitcellhull==ConvexHull[Join[unitcellvert,{til[[6,i]]}]],AppendTo[indices,i]],
		{i,Range[n]}];

		Print[Length[indices]," tiles in unit cell"];

		tiles = til[[4,indices]];
		tilpoints = Union[Flatten[tiles,1]];
		
	{indices,tiles,tilpoints,unitcellvert,mapvert,basisvectors}
]


DualizeUnitCell[unitcell_,starvectors_,modvectors_,a_:0,b_:0,
			ptsize_:0.05,
			linewidth_:1/200]:=
	Module[{SM,indices,tiles,tilpoints,unitcellvert,mapvert,basisvectors,centerpolygons,GetNewBasis,b1,b2,R,Rinv,com,A,Ainv,translate,bonds},

		{indices,tiles,tilpoints,unitcellvert,mapvert,basisvectors}=unitcell;

		SM={{1,a},{b,1}};
		
		centerpolygons=N[Map[add[#,modvectors]&,tiles,{2}]];

		centerpolygons=Map[SM.#&,centerpolygons,{2}];

		mapvert=SM.#&/@mapvert;

		GetNewBasis[{basis1_,basis2_}]:=
			Module[{t,lattice,latticepts,norms,newa1,counter,found,curpt,newa2},
				t=Range[-4,4];
				lattice=Table[{t1,t2,t1*basis1+t2*basis2},{t1,t},{t2,t}];
				latticepts=Flatten[lattice,1];
				norms=Norm/@latticepts[[All,3]];
				newa1=latticepts[[Position[norms,RankedMin[norms,2]][[-1]]]][[1]];		
				counter=3;
				found=False;
				While[!found,
					curpt=latticepts[[Position[norms,RankedMin[norms,counter]][[-1]]]][[1]];
					If[Chop[N[Norm[Cross[Join[newa1[[3]],{0}],Join[curpt[[3]],{0}]]]]]!=0,
					newa2=curpt;
					found=True;,
					counter++;]
				];
				newa2=curpt;
			{newa1[[3]],newa2[[3]]}
		];

		basisvectors=GetNewBasis[SM.#&/@basisvectors];

		(*define b1,b2 so that b1 has larger y-value*)
		If[basisvectors[[1,2]]>basisvectors[[2,2]],
			b1=basisvectors[[1]]; b2=basisvectors[[2]],
			b1=basisvectors[[2]]; b2=basisvectors[[1]]]; 

		mapvert={mapvert[[4]]+b2,mapvert[[4]]+b1+b2,mapvert[[4]]+b1,mapvert[[4]]};

		(*rotate everything so that b1 points vertically*)
		R=RotationMatrix[{Normalize[b1],{0,1}}];
		Rinv=Inverse[R];
		basisvectors=Chop[Map[R.#&,{b1,b2}]];
		mapvert=Map[R.#&,mapvert];
				
		(*translate by looking at skewed coordinates, taking mod to save only those in unit
		cell, then adding translation vector to mapvert lower left*)

		com=mapvert[[4]];
		A=Transpose[{basisvectors[[1]],basisvectors[[2]]}];
		Ainv=Inverse[A];

		centerpolygons=Map[R.#&,centerpolygons,{2}];
		centerpolygons=Map[#-com&,centerpolygons,{2}];
		centerpolygons=Map[Ainv.#&,centerpolygons,{2}];
		translate=Map[-Total[#]/4+Mod[Total[#]/4,1]&,centerpolygons];
		centerpolygons=centerpolygons+Transpose[{translate,translate,translate,translate}];
		centerpolygons=Map[A.#&,centerpolygons,{2}];
		centerpolygons=Map[com+#&,centerpolygons,{2}];
		centerpolygons=Map[#-mapvert[[4]]&,centerpolygons,{2}];
		mapvert=Map[#-mapvert[[4]]&,mapvert];

		bonds=Flatten[Map[{{#[[1]],#[[2]]},{#[[2]],#[[3]]},{#[[3]],#[[4]]},{#[[4]],#[[1]]}}&,Chop[centerpolygons]],1];
		bonds=DeleteDuplicates[bonds];

	Chop[{centerpolygons,mapvert,basisvectors,bonds}]
]


GetInfoGraphics[acc_,uc_,maxradius_,newangles_]:=
	Module[{CurAllAngles=newangles,Angles,numpoints,infogfx},
		AppendTo[CurAllAngles,360-Total[CurAllAngles]];
		Angles=StringJoin["\!\(\*SubscriptBox[\(\[Theta]\), \(ideal\)]\) = ", ToString[CurAllAngles]];
		numpoints=Length[uc[[1]]];

		infogfx=Show[Graphics[Text[Style[StringJoin[
			Angles,"\n",
			"Accuracy = ",ToString[acc],"\n",
			"Num. of points = ",ToString[numpoints],"\n",
			"Min. distance / 2 (max radius) = ",ToString[maxradius],"\n"]
		,FontFamily->"Times",FontSize->30]],PlotRange->{{-1,1},{-1,1}}]];
	infogfx
]


GetKLattice[basisvector1_,basisvector2_]:=
	Module[{R90,b1,b2,CalculateTwoHexFaces,G,M0,M1,M2,K0,K1,K2,l1,v1,l2,v2,l3,v3,Kpoints},

		R90=RotationMatrix[Pi/2];
		b1=2 Pi R90.basisvector2 / (basisvector1.(R90.basisvector2));
		b2=2 Pi R90.basisvector1 / (basisvector2.(R90.basisvector1));

		(*calculates the norm and direction for the K point, given vectors touching centers of neighboring edges (in coordinates of basis vectors, which are also given)*)
		CalculateTwoHexFaces[v1_,v2_,b1_,b2_]:=
			Module[{n1,n2,a1,a2,alpha,roots,roots2,n,x,theta,a,b,c},
				a1=v1.{b1,b2}; a2=v2.{b1,b2};
				n1=Norm[a1]; n2=Norm[a2];
				alpha=ArcCos[Normalize[a1].Normalize[a2]];
				Clear[x];
				roots=Solve[Cos[alpha]n2/x+Sin[alpha]Sqrt[1-(n2/x)^2]==n1/x,x];
				n=x/.roots[[1]];
				roots2=Solve[1==x^2n1^2+n2^2((n1^2/n-x n1^2)/Cos[alpha]/n1/n2)^2+2x(n1^2/n-x n1^2),x];
				a=x/.roots2[[1]];
				b=(n1^2/n-a n1^2)/Cos[alpha]/n1/n2;
				c=a v1 + b v2;
			{n,c}
		];

		If[ArcCos[Normalize[b1].Normalize[b2]]>Pi/2,
			M0={0.5,0};
			M1={0.5,0.5};
			M2={0,0.5};,
			
			M0={0.5,0};
			M1={0,0.5};
			M2={-0.5,0.5};];
		G={0,0};

		{l1,v1}=CalculateTwoHexFaces[M0,M1,b1,b2];
		{l2,v2}=CalculateTwoHexFaces[M1,M2,b1,b2];
		{l3,v3}=CalculateTwoHexFaces[M2,-M0,b1,b2];
		K0=l1 v1;
		K1=l2 v2;
		K2=l3 v3;

		Kpoints={G,M0,K0,M1,K1,M2,K2};

	{{b1,b2},Kpoints}
]


PlotBZ[Kcell_]:=
	Module[{Kbasis,Kpoints,recipgfx,t,klattice,klatticegfx,a,scale,xmin,ymin,xmax,ymax,pklsize,pkl},
		{Kbasis,Kpoints}=Kcell;

		a=Join[Kbasis,Kbasis*-1];
		scale=0.1;
		{xmin,ymin}=(1+scale)*Min/@Transpose[a];
		{xmax,ymax}=(1+scale)*Max/@Transpose[a];
		pklsize={{xmin,xmax},{ymin,ymax}};

		recipgfx=Graphics[{Red,Map[Arrow[{{0,0},#}]&,Kbasis]}];

		t=Range[-2,2];
		klattice=Flatten[Table[t1*Kbasis[[1]]+t2*Kbasis[[2]],{t1,t},{t2,t}],1];

		klatticegfx=Graphics[Point/@klattice];

		pkl=Show[{VoronoiMesh[klattice],klatticegfx,recipgfx,Graphics[Arrow[{{0,0},Total[Kbasis]/2}]],Graphics[{PointSize[Large],Point/@Map[#.Kbasis&,Kpoints]}]},
			PlotRange->pklsize];
	pkl
]


TrimUnitCell3[unitcell_,starvectors_,modvectors_,a_:0,b_:0]:=
	Module[{SM,indices,tiles,tilpoints,unitcellvert,mapvert,basisvectors,sv,GetArea,GetNewBasis,b1,b2,R,Rinv,basisgfx,
		pts,com,A,Ainv,tol,
		interiorindices={},
		translatedpoints={}
		},

		SM={{1,a},{b,1}};

		{indices,tiles,tilpoints,unitcellvert,mapvert,basisvectors}=unitcell;
		{mapvert,basisvectors}=N/@{mapvert,basisvectors};

		mapvert=SM.#&/@mapvert;
		

		GetArea[v1_,v2_]:=Norm[Cross[Join[v1,{0}],Join[v2,{0}]]];

		GetNewBasis[{basis1_,basis2_}]:=
			Module[{t,lattice,latticepts,norms,newa1,counter,found,curpt,newa2},
				t=Range[-4,4];
				lattice=Table[{t1,t2,t1*basis1+t2*basis2},{t1,t},{t2,t}];
				latticepts=Flatten[lattice,1];
				norms=Norm/@latticepts[[All,3]];
				newa1=latticepts[[Position[norms,RankedMin[norms,2]][[-1]]]][[1]];		
				counter=3;
				found=False;
				While[!found,
					curpt=latticepts[[Position[norms,RankedMin[norms,counter]][[-1]]]][[1]];
					If[Chop[N[Norm[Cross[Join[newa1[[3]],{0}],Join[curpt[[3]],{0}]]]]]!=0,
					newa2=curpt;
					found=True;,
					counter++;]
				];
				newa2=curpt;
			{newa1[[3]],newa2[[3]]}
		];

		(*need to choose different basis*)
		basisvectors=GetNewBasis[SM.#&/@basisvectors];
		
		(*define b1,b2 so that b1 has larger y-value*)
		If[basisvectors[[1,2]]>basisvectors[[2,2]],
			b1=basisvectors[[1]]; b2=basisvectors[[2]],
			b1=basisvectors[[2]]; b2=basisvectors[[1]]]; 

		mapvert={mapvert[[4]]+b2,mapvert[[4]]+b1+b2,mapvert[[4]]+b1,mapvert[[4]]};

		(*rotate everything so that b1 points vertically*)
		R=RotationMatrix[{Normalize[b1],{0,1}}];
		Rinv=Inverse[R];
		basisvectors=Chop[Map[R.#&,{b1,b2}]];
		mapvert=Map[R.#&,mapvert];

		basisgfx={Line[{mapvert[[4]],mapvert[[3]]}],
					Line[{mapvert[[4]],mapvert[[1]]}],
					Line[{mapvert[[1]],mapvert[[2]]}],
					Line[{mapvert[[3]],mapvert[[2]]}]};

		(*sv = R.#&/@starvectors[[2]];
		pts = Map[add[#,sv]&,tilpoints];
		pts = SM.#&/@pts;*)

		
		pts = Map[add[#,modvectors]&,tilpoints];
		pts = SM.#&/@pts;
		pts = R.#&/@pts;
		
				
		(*translate by looking at skewed coordinates, taking mod to save only those in unit
		cell, then adding translation vector to mapvert lower left*)

		com=mapvert[[4]];
		pts=#-com&/@pts;
		A=Transpose[{basisvectors[[1]],basisvectors[[2]]}];
		Ainv=Inverse[A];
		pts=Ainv.#&/@pts;
		pts=Mod[#,1]&/@pts;
		pts=A.#&/@pts;
		pts=#+com&/@pts;

		(*Print[Show[Graphics[basisgfx],Graphics[Point/@pts]]];*)
		tol=0.005;
		pts = DeleteCases[pts,_?(Norm[#-mapvert[[1]]]< tol ||
								Norm[#-mapvert[[2]]]< tol ||
								Norm[#-mapvert[[4]]]< tol&),{1},Heads->False];

		AppendTo[pts,mapvert[[3]]];
		pts = DeleteDuplicates[pts,Norm[#1-#2]<tol &];

		pts=#-mapvert[[4]]&/@pts;
		mapvert=Map[#-mapvert[[4]]&,mapvert];
	
		Print[Length[pts]," points in trimmed, translated unit cell"];
		
		{pts,mapvert,basisvectors}
]


PlotDualUnitCell[dualcell_,
			ptsize_:0.05,
			linewidth_:1/200]:=
	Module[{centerpolygons,mapvert,basisvectors,bonds,basisgfx,centerpolygonsgfx,bondsgfx,
		NN1,NN2,polyNN1,polyNN2,polyNN1gfx,polyNN2gfx,bondsNN1,bondsNN1gfx,bondsNN2,bondsNN2gfx,
		xmin,ymin,xmax,ymax,xlen,ylen,pxmin,pxmax,pymin,pymax,scale,pucsize,puc},
		
		{centerpolygons,mapvert,basisvectors,bonds} = dualcell[[{1,2,3,4}]];

		basisgfx={Line[{mapvert[[4]],mapvert[[3]]}],
					Line[{mapvert[[4]],mapvert[[1]]}],
					Line[{mapvert[[1]],mapvert[[2]]}],
					Line[{mapvert[[3]],mapvert[[2]]}]};

		centerpolygonsgfx = {Black,Opacity[0.5],Polygon/@centerpolygons};

		NN1=Table[Apply[Plus,x*basisvectors],{x,{{1,0},{-1,0},{0,1},{0,-1}}}];
		NN2=Table[Apply[Plus,x*basisvectors],{x,{{1,1},{1,-1},{-1,1},{-1,-1}}}];

		polyNN1 = Table[Table[Map[#+x&,y],{y,centerpolygons}],{x,NN1}];
		polyNN2 = Table[Table[Map[#+x&,y],{y,centerpolygons}],{x,NN2}];

		polyNN1gfx = Map[{Blue,Opacity[0.5],Polygon[#]}&,polyNN1,{2}];
		polyNN2gfx = Map[{Orange,Opacity[0.5],Polygon[#]}&,polyNN2,{2}];

		bondsNN1=Table[Table[Map[#+x&,y],{y,bonds}],{x,NN1}];
		bondsNN2=Table[Table[Map[#+x&,y],{y,bonds}],{x,NN2}];

		bondsgfx=Line/@bonds;
		bondsNN1gfx=Map[{Black,Opacity[0.2],Line[#]}&,bondsNN1];
		bondsNN2gfx=Map[{Black,Opacity[0.2],Line[#]}&,bondsNN2];

		{xmin,ymin}=Min/@Transpose[mapvert];
		{xmax,ymax}=Max/@Transpose[mapvert];
		xlen=xmax-xmin;
		ylen=ymax-ymin;

		scale=6;
		pxmin=xmin-xlen/scale;
		pxmax=xmax+xlen/scale;
		pymin=ymin-ylen/scale;
		pymax=ymax+ylen/scale;
		pucsize={{pxmin,pxmax},{pymin,pymax}};

		puc=Graphics[{basisgfx,{PointSize[0.02],Map[Point,mapvert]},
						centerpolygonsgfx,polyNN1gfx,polyNN2gfx,bondsNN1gfx,bondsNN2gfx,bondsgfx},PlotRange->pucsize];
	puc
]


GetIntersections[til_,v_]:=
Module[{vi,rows,intersecpts={}},
		Do[
			rows=Join[Position[til[[5]],{vi[[1]],vi[[2]]}],Position[til[[5]],{vi[[2]],vi[[1]]}]][[All,1]];
			intersecpts=Join[intersecpts,Part[til[[6]],rows]];,
			{vi,v}
		];
	intersecpts
]


GetIntersectionTiles[kmin_,
					kmax_,
					F1_,F0_,
					v_,
					gamma_:"random"
					gammaseed_:1]:=
Module[{r,rorth,rmag,nj,kmn,kmx,gam,cut,
		i1,i2,origin,pre,b1,b2,
		intersectpts,tiles,tillines},

		{r,rorth,rmag} = StarVectors[F1,F0];
		Print["star vector magnitudes:", N[rmag]];
		SeedRandom[gammaseed];

		nj         = 5;
		kmn        = Which[And[Head[kmin]===List,
                                 Length[kmin]===nj,
                                 Union[Map[IntegerQ,kmin]]==={True}],
                             kmin,
                             IntegerQ[kmin],
                             Table[kmin,{nj}],
                             True,
                             Print["Error: invalid argument kmin"];
                             Return[]];
          kmx        = Which[And[Head[kmax]===List,
                                 Length[kmax]===nj,
                                 Union[Map[IntegerQ,kmax]]==={True}],
                             kmax,
                             IntegerQ[kmax],
                             Table[kmax,{nj}],
                             True,
                             Print["Error: invalid argument kmax"];
                             Return[]];
          gam        = Which[gamma==="random",
                         Print["note: translation vector gamma ",
                               "is chosen randomly!"];
                         Table[Random[Real],{nj}],
						gamma==="rand2",
						Print["note: translation vector gamma ",
							  "is chosen for first two components"];
						StarPhases[F1,F0],
                         And[Head[gamma]===List,
                             Length[gamma]==2,
                             First[gamma]==="random",
                             NumberQ[Last[gamma]]],
                         Print["note: translation vector gamma is chosen"];
                         Print["      randomly (sum of components is ",
                               gamma[[2]],")"];
                         Apply[Append[#,gamma[[2]]-Apply[Plus,#]]&,
                               {Table[Random[Real,{-0.5,0.5}],{nj-1}]}],
                         And[Head[gamma]===List,
                             Length[gamma]===nj,
                             Union[Map[NumberQ,gamma]]==={True}],
                         gamma,
                         NumberQ[gamma],
                         Table[gamma,{nj}],
                         True,
                         Print["Error: invalid argument gamma"];
                         Return[]];


		i1=v[[1]];
		i2=v[[2]];
		origin=1/Dot[rorth[[i1]],r[[i2]]]*(rmag[[i1]]gam[[i1]]rorth[[i2]]-
										rmag[[i2]]gam[[i2]]rorth[[i1]]);
		pre=1/Sqrt[rmag[[i1]]^2rmag[[i2]]^2-Dot[r[[i1]],r[[i2]]]^2];
		b1=pre*rmag[[i2]]^2 rorth[[i1]];
		b2=pre*rmag[[i1]]^2 rorth[[i2]];

		intersectpts = Flatten[Table[origin+k1*b1+k2*b2,
									{k1,kmn[[i1]],kmx[[i1]]},{k2,kmn[[i2]],kmx[[i2]]}]
								,1];
		
		tiles = Table[{x,x+b1,x+b1+b2,x+b2},{x,intersectpts}];
		tillines  = Flatten[Map[Partition[#,2,1]&,
                                  Map[Append[#,#[[1]]]&,
                                      tiles]],1];
		tillines  = Union[Map[Sort,tillines]];
          
	{intersectpts,tiles,tillines}
]

PlotIntersectionTiles[kmin_,
					kmax_,
					F1_,F0_,
					v_,
					gamma_:"random",
					theta_:2Pi/5,
					showpolygons_:True,
					showlines_:True,
					showpoints_:False,
					ptsize_:1/100,
					linewidth_:1/200,
					plotrange_:6]:=
Module[{g,vi,ti,points={},lines={},polygons={},
		plotlines,plotpoints,plotpolygons},

	g = PlotGrid[kmin,kmax,F1,F0,gamma];
	
	Do[
		ti = GetIntersectionTiles[kmin,kmax,F1,F0,vi,gamma];
		AppendTo[points,Table[Join[{Green,Opacity[0.5],PointSize[ptsize]},Map[Point,{x}]],{x,ti[[1,All]]}]];
		AppendTo[polygons,Table[Join[{Green,Opacity[0.15]},Map[Polygon,{x}]],{x,ti[[2,All]]}]];
		AppendTo[lines,Table[Join[{Thickness[linewidth],Green},Map[Line,{x}]],{x,ti[[3,All]]}]];,
		{vi,v}];
	plotpoints = If[showpoints,Flatten[points],{}];
	plotlines = If[showlines,Flatten[lines],{}];
	plotpolygons = If[showpolygons,Flatten[polygons],{}];


	Show[Graphics[Join[plotpoints,plotpolygons,plotlines]],
			AspectRatio -> Automatic,
			PlotRange->plotrange{{-1,1},{-1,1}}]
]


TrimUnitCell[unitcell_,F1_,F0_,theta_:2Pi/5]:=
Module[{indices,tiles,tilpoints,unitcellvert,mapvert,basisvectors,b1,b2,
		R,r,R2,R3,
		pts,unitcellhull,n,tol,
		leftind,topind,rightind,bottomind,yind,xind,cind,
		m,b,corners,
		interiorindices={},
		translatedpoints={},
		MoveToInterior,
		basisvectorsgfx},

		{indices,tiles,tilpoints,unitcellvert,mapvert,basisvectors}=unitcell;

		mapvert=N[mapvert];
		basisvectors=N[basisvectors];

		If[basisvectors[[1,2]]>basisvectors[[2,2]],
			b1=basisvectors[[1]]; b2=basisvectors[[2]],
			b1=basisvectors[[2]]; b2=basisvectors[[1]]]; (*b1 has larger y-value*)
		R=RotationMatrix[{Normalize[b1],{0,1}}];


		basisvectors=Chop[Map[R.#&,basisvectors]];
		mapvert=Map[R.#&,mapvert];
		r = StarVectors[F1,F0,theta];
		r[[1]]=Table[R.x,{x,Map[Normalize,r[[1]]]}];

		pts = Map[add2[#,r]&,tilpoints];
		
		n = Length[tilpoints];

		(*
		Print[mapvert];
		Print[basisvectors];
		*)

		(*moves all boundary points to left and top*)
		(*assumes mapvert[[1]] is on right side*)
		(*tol is tolerance for how close point is to edge*)
		tol=0.05;

		(*translate left points*)
		leftind=Flatten[Position[pts,_?(#[[1]]<mapvert[[3,1]]-tol &),{1},Heads->False]];
		Do[pts[[i]]=pts[[i]]+basisvectors[[2]],
			{i,leftind}];

		(*translate right boundary points*)
		(*rightind=Flatten[Position[pts,
								_?(Abs[#[[1]]-mapvert[[1,1]]]< tol &),{1},Heads->False]];*)
		rightind=Flatten[Position[pts,_?(#[[1]]>mapvert[[1,1]]- tol &),{1},Heads->False]];
		Do[pts[[i]]=pts[[i]]-basisvectors[[2]],
			{i,rightind}];

		(*translate top points*)
		m = (mapvert[[1,2]]-mapvert[[4,2]])/(mapvert[[1,1]]-mapvert[[4,1]]);
		b = mapvert[[2,2]] - m * mapvert[[2,1]];
		topind=Flatten[Position[pts,_?(#[[2]]> m*#[[1]]+b &),{1},Heads->False]];
		Do[pts[[i]]=pts[[i]]-basisvectors[[1]],
			{i,topind}];

		(*translate bottom boundary points*)
		(*line defining bottom edge, y = mx + b*)
		b = mapvert[[1,2]] - m * mapvert[[1,1]];

		(*bottomind=Flatten[Position[pts,
								_?(Abs[#[[2]]-(m*#[[1]]+b)]< tol &),{1},Heads->False]];*)
		bottomind=Flatten[Position[pts,_?(#[[2]]<m*#[[1]]+b+ tol &),{1},Heads->False]];
		Do[pts[[i]]=pts[[i]]+basisvectors[[1]],
			{i,bottomind}];

		(*leave only top left corner, assuming specified by mapvert[[3]]*)
		pts = DeleteCases[pts,_?(Norm[#-mapvert[[1]]]< tol ||
								Norm[#-mapvert[[2]]]< tol ||
								Norm[#-mapvert[[4]]]< tol&),{1},Heads->False];

		pts = DeleteDuplicates[pts,Norm[#1-#2]<tol &];

		(*ABOVE IS ALL FROM PREVIOUS CODE... NOW DO EXTRA CHECKS FOR ROTATED BASIS, i.e. a basis vector not perpendicular with y axis*)

(*
		yind=Flatten[Position[pts,_?(#[[2]]>mapvert[[3,2]]- tol &),{1},Heads->False]];
		Do[pts[[i]]=pts[[i]]+basisvectors[[2]],
			{i,yind}];

		xind=Flatten[Position[pts,_?(#[[1]]<mapvert[[4,1]]- tol &),{1},Heads->False]];
		Do[pts[[i]]=pts[[i]]+basisvectors[[2]],
			{i,xind}];

		cind=Flatten[Position[pts,_?(Cross[Append[#-mapvert[[4]],0],Append[basisvectors[[1]],0]][[3]]<0 &),{1},Heads->False]];
		Do[pts[[i]]=pts[[i]]+basisvectors[[2]],
			{i,cind}];
*)
(*leave only top left corner, assuming specified by mapvert[[3]]*)
		pts = DeleteCases[pts,_?(Norm[#-mapvert[[1]]]< tol ||
								Norm[#-mapvert[[2]]]< tol ||
								Norm[#-mapvert[[4]]]< tol&),{1},Heads->False];


		AppendTo[pts,mapvert[[3]]];
		pts = DeleteDuplicates[pts,Norm[#1-#2]<tol &];

	Print[Length[pts]," points in unit cell"];
	pts
]

TrimUnitCell2[unitcell_,angles_,theta_,indexlist_]:=
Module[{indices,tiles,tilpoints,unitcellvert,mapvert,basisvectors,b1,b2,
		R,r,R2,R3,
		pts,unitcellhull,n,tol,
		leftind,topind,rightind,bottomind,yind,xind,cind,
		m,b,corners,t,
		interiorindices={},
		translatedpoints={},
		MoveToInterior,
		basisvectorsgfx},

		{indices,tiles,tilpoints,unitcellvert,mapvert,basisvectors}=unitcell;

		mapvert=N[mapvert];
		basisvectors=N[basisvectors];

		If[basisvectors[[1,2]]>basisvectors[[2,2]],
			b1=basisvectors[[1]]; b2=basisvectors[[2]],
			b1=basisvectors[[2]]; b2=basisvectors[[1]]]; (*b1 has larger y-value*)
		R=RotationMatrix[{Normalize[b1],{0,1}}];

		r = StarVectors2[angles,theta,indexlist];
		r[[1]]=Normalize/@r[[1]];

		pts = Map[add2[#,r]&,tilpoints];

		t=mapvert[[4]];
		pts=#-t&/@pts;
		mapvert=#-t&/@mapvert;

		basisvectors=Chop[Map[R.#&,basisvectors]];
		mapvert=Map[R.#&,mapvert];
		pts=R.#&/@pts;
		
		n = Length[tilpoints];

		(*
		Print[mapvert];
		Print[basisvectors];
		*)

		(*moves all boundary points to left and top*)
		(*assumes mapvert[[1]] is on right side*)
		(*tol is tolerance for how close point is to edge*)
		tol=0.005;

		(*translate left points*)
		leftind=Flatten[Position[pts,_?(#[[1]]<mapvert[[3,1]]-tol &),{1},Heads->False]];
		Do[pts[[i]]=pts[[i]]+basisvectors[[2]],
			{i,leftind}];

		(*translate right boundary points*)
		(*rightind=Flatten[Position[pts,
								_?(Abs[#[[1]]-mapvert[[1,1]]]< tol &),{1},Heads->False]];*)
		rightind=Flatten[Position[pts,_?(#[[1]]>mapvert[[1,1]]- tol &),{1},Heads->False]];
		Do[pts[[i]]=pts[[i]]-basisvectors[[2]],
			{i,rightind}];

		(*translate top points*)
		m = (mapvert[[1,2]]-mapvert[[4,2]])/(mapvert[[1,1]]-mapvert[[4,1]]);
		b = mapvert[[2,2]] - m * mapvert[[2,1]];
		topind=Flatten[Position[pts,_?(#[[2]]> m*#[[1]]+b &),{1},Heads->False]];
		Do[pts[[i]]=pts[[i]]-basisvectors[[1]],
			{i,topind}];

		(*translate bottom boundary points*)
		(*line defining bottom edge, y = mx + b*)
		b = mapvert[[1,2]] - m * mapvert[[1,1]];

		(*bottomind=Flatten[Position[pts,
								_?(Abs[#[[2]]-(m*#[[1]]+b)]< tol &),{1},Heads->False]];*)
		bottomind=Flatten[Position[pts,_?(#[[2]]<m*#[[1]]+b+ tol &),{1},Heads->False]];
		Do[pts[[i]]=pts[[i]]+basisvectors[[1]],
			{i,bottomind}];

		(*leave only top left corner, assuming specified by mapvert[[3]]*)
		pts = DeleteCases[pts,_?(Norm[#-mapvert[[1]]]< tol ||
								Norm[#-mapvert[[2]]]< tol ||
								Norm[#-mapvert[[4]]]< tol&),{1},Heads->False];

		pts = DeleteDuplicates[pts,Norm[#1-#2]<tol &];

		(*ABOVE IS ALL FROM PREVIOUS CODE... NOW DO EXTRA CHECKS FOR ROTATED BASIS, i.e. a basis vector not perpendicular with y axis*)

(*
		yind=Flatten[Position[pts,_?(#[[2]]>mapvert[[3,2]]- tol &),{1},Heads->False]];
		Do[pts[[i]]=pts[[i]]+basisvectors[[2]],
			{i,yind}];

		xind=Flatten[Position[pts,_?(#[[1]]<mapvert[[4,1]]- tol &),{1},Heads->False]];
		Do[pts[[i]]=pts[[i]]+basisvectors[[2]],
			{i,xind}];

		cind=Flatten[Position[pts,_?(Cross[Append[#-mapvert[[4]],0],Append[basisvectors[[1]],0]][[3]]<0 &),{1},Heads->False]];
		Do[pts[[i]]=pts[[i]]+basisvectors[[2]],
			{i,cind}];
*)
(*leave only top left corner, assuming specified by mapvert[[3]]*)
		pts = DeleteCases[pts,_?(Norm[#-mapvert[[1]]]< tol ||
								Norm[#-mapvert[[2]]]< tol ||
								Norm[#-mapvert[[4]]]< tol&),{1},Heads->False];


		AppendTo[pts,mapvert[[3]]];
		pts = DeleteDuplicates[pts,Norm[#1-#2]<tol &];

	Print[Length[pts]," points in unit cell"];
	pts
]



PlotSingular[til_,
			unitcell_,
			kmin_,
			kmax_,
			gamma_,
			F1_,
			F0_,
			tol_,
			plotrange_:6]:=
Module[{indices,tiles,tilpoints,mapvert,basisvectors,unitcellvert,unitcellbdry,tip},
		
	{indices,tiles,tilpoints,unitcellvert,mapvert,basisvectors}=unitcell;
	
	unitcellbdry = {Green,Opacity[0.5],Polygon[unitcellvert]};

	tip=Tally[til[[6]],Norm[#1-#2]<tol&];
	Show[Graphics[
			{unitcellbdry,PlotGrid[kmin,kmax,F1,F0,gamma],
			{Cyan,PointSize[Large],Map[Point,til[[6,indices]]]},
			{Blue,PointSize[Large],Map[Point,tip[[Flatten[Position[tip[[All,2]],5]],1]]]},
			{Red,PointSize[Large],Map[Point,tip[[Flatten[Position[tip[[All,2]],4]],1]]]},
			{Green,PointSize[Large],Map[Point,tip[[Flatten[Position[tip[[All,2]],3]],1]]]},
			{Orange,PointSize[Large],Map[Point,tip[[Flatten[Position[tip[[All,2]],2]],1]]]},
			{Black,PointSize[0.01],Map[Point,tip[[Flatten[Position[tip[[All,2]],1]],1]]]}}],
		AspectRatio->Automatic,
		PlotRange->plotrange{{-1,1},{-1,1}}
	]
]

PlotColors[til_,unitcell_,F1_,F0_,
			groups_:{1,2,3,4},
			showpolygons_:True,
			showlines_:True,
			showpoints_:True,
			linewidth_:1/500,
			ptsize_:1/100,
			len_:20] :=
Module[{r,rorth,rmag,g,
		indices,tiles,tilpoints,unitcellvert,mapvert,basisvectors,
		g12,g13,g14,g15,g23,g24,g25,g34,g35,g45,
		points1,points2,points3,points4,lines1,lines2,lines3,lines4,polygons1,polygons2,polygons3,polygons4,
		unitpolygons,
		points,lines,polygons,
		plotlines,plotpoints,plotpolygons},

		{r,rorth,rmag} = StarVectors[F1,F0];

		g = Table[{r[[j]] - len*rorth[[j]],
                    r[[j]] + len*rorth[[j]]},
                 {j,2}];
		g = Map[Line,g];

		{indices,tiles,tilpoints,unitcellvert,mapvert,basisvectors}=unitcell;

		g12=GetSingleTile[til,1,2];
		g13=GetSingleTile[til,1,3];
		g14=GetSingleTile[til,1,4];
		g15=GetSingleTile[til,1,5];
		g23=GetSingleTile[til,2,3];
		g24=GetSingleTile[til,2,4];
		g25=GetSingleTile[til,2,5];
		g34=GetSingleTile[til,3,4];
		g35=GetSingleTile[til,3,5];
		g45=GetSingleTile[til,4,5];

		points1 = Table[Join[{PointSize[ptsize],Green,Opacity[0.8]},Map[Point,GetMainVertex[x]]],{x,{g12,g13,g14,g15}}];
		points2 = Table[Join[{PointSize[ptsize],Blue,Opacity[0.8]},Map[Point,GetMainVertex[x]]],{x,{g23,g24,g25}}];
		points3 = Table[Join[{PointSize[ptsize],Red,Opacity[0.8]},Map[Point,GetMainVertex[x]]],{x,{g34,g35}}];
		points4 = Table[Join[{PointSize[ptsize],Black,Opacity[0.8]},Map[Point,GetMainVertex[x]]],{x,{g45}}];

		points={points1,points2,points3,points4};

		lines1 = Table[Join[{Thickness[linewidth],Green},Map[Line,GetLines[x]]],{x,{g12,g13,g14,g15}}];
		lines2 = Table[Join[{Thickness[linewidth],Blue},Map[Line,GetLines[x]]],{x,{g23,g24,g25}}];
		lines3 = Table[Join[{Thickness[linewidth],Red},Map[Line,GetLines[x]]],{x,{g34,g35}}];
		lines4 = Table[Join[{Thickness[linewidth],Black},Map[Line,GetLines[x]]],{x,{g45}}];

		lines={lines1,lines2,lines3,lines4};

		polygons1 = Table[Join[{Green,Opacity[0.5]},Map[Polygon,GetTiles[x]]],{x,{g12,g13,g14,g15}}];
		polygons2 = Table[Join[{Blue,Opacity[0.5]},Map[Polygon,GetTiles[x]]],{x,{g23,g24,g25}}];
		polygons3 = Table[Join[{Red,Opacity[0.5]},Map[Polygon,GetTiles[x]]],{x,{g34,g35}}];
		polygons4 = Table[Join[{Black,Opacity[0.5]},Map[Polygon,GetTiles[x]]],{x,{g45}}];

		unitpolygons = Map[{Purple,Opacity[0.5],#}&,Map[Polygon,Map[add,tiles,{2}]]];

		polygons={polygons1,polygons2,polygons3,polygons4};

		plotlines = If[showlines,Flatten[lines[[groups]]],{}];
		plotpoints = If[showpoints,Flatten[points[[groups]]],{}];
		plotpolygons = If[showpolygons,Flatten[polygons[[groups]]],{}];
		Show[Graphics[Join[g,plotlines,plotpolygons,plotpoints,unitpolygons]],
			AspectRatio -> Automatic]
];

PlotColors2[til_,F1_,F0_,
			groups_:{1,2,3,4},
			showpolygons_:True,
			showlines_:True,
			showpoints_:True,
			linewidth_:1/500,
			ptsize_:1/100,
			len_:20] :=
Module[{r,rorth,rmag,g,
		g12,g13,g14,g15,g23,g24,g25,g34,g35,g45,
		points1,points2,points3,points4,lines1,lines2,lines3,lines4,polygons1,polygons2,polygons3,polygons4,
		points,lines,polygons,
		plotlines,plotpoints,plotpolygons},

		{r,rorth,rmag} = StarVectors[F1,F0];

		g = Table[{r[[j]] - len*rorth[[j]],
                    r[[j]] + len*rorth[[j]]},
                 {j,2}];
		g = Map[Line,g];
		
		g12=GetSingleTile[til,1,2];
		g13=GetSingleTile[til,1,3];
		g14=GetSingleTile[til,1,4];
		g15=GetSingleTile[til,1,5];
		g23=GetSingleTile[til,2,3];
		g24=GetSingleTile[til,2,4];
		g25=GetSingleTile[til,2,5];
		g34=GetSingleTile[til,3,4];
		g35=GetSingleTile[til,3,5];
		g45=GetSingleTile[til,4,5];

		points1 = Table[Join[{PointSize[ptsize],Green,Opacity[0.8]},Map[Point,GetMainVertex2[x,F1,F0]]],{x,{g12,g13,g14,g15}}];
		points2 = Table[Join[{PointSize[ptsize],Blue,Opacity[0.8]},Map[Point,GetMainVertex2[x,F1,F0]]],{x,{g23,g24,g25}}];
		points3 = Table[Join[{PointSize[ptsize],Red,Opacity[0.8]},Map[Point,GetMainVertex2[x,F1,F0]]],{x,{g34,g35}}];
		points4 = Table[Join[{PointSize[ptsize],Black,Opacity[0.8]},Map[Point,GetMainVertex2[x,F1,F0]]],{x,{g45}}];

		points={points1,points2,points3,points4};

		lines1 = Table[Join[{Thickness[linewidth],Green},Map[Line,GetLines2[x,F1,F0]]],{x,{g12,g13,g14,g15}}];
		lines2 = Table[Join[{Thickness[linewidth],Blue},Map[Line,GetLines2[x,F1,F0]]],{x,{g23,g24,g25}}];
		lines3 = Table[Join[{Thickness[linewidth],Red},Map[Line,GetLines2[x,F1,F0]]],{x,{g34,g35}}];
		lines4 = Table[Join[{Thickness[linewidth],Black},Map[Line,GetLines2[x,F1,F0]]],{x,{g45}}];

		lines={lines1,lines2,lines3,lines4};

		polygons1 = Table[Join[{Green,Opacity[0.5]},Map[Polygon,GetTiles2[x,F1,F0]]],{x,{g12,g13,g14,g15}}];
		polygons2 = Table[Join[{Blue,Opacity[0.5]},Map[Polygon,GetTiles2[x,F1,F0]]],{x,{g23,g24,g25}}];
		polygons3 = Table[Join[{Red,Opacity[0.5]},Map[Polygon,GetTiles2[x,F1,F0]]],{x,{g34,g35}}];
		polygons4 = Table[Join[{Black,Opacity[0.5]},Map[Polygon,GetTiles2[x,F1,F0]]],{x,{g45}}];

		polygons={polygons1,polygons2,polygons3,polygons4};

		plotlines = If[showlines,Flatten[lines[[groups]]],{}];
		plotpoints = If[showpoints,Flatten[points[[groups]]],{}];
		plotpolygons = If[showpolygons,Flatten[polygons[[groups]]],{}];

		Show[Graphics[Join[g,plotlines,plotpolygons,plotpoints]],
			AspectRatio -> Automatic]
];


GetPoints[til_]:= Map[add,til[[2]]]

GetTiles[til_]:=Map[add,til[[4]],{2}]
GetLines[til_]:=Map[add,til[[3]],{2}]
GetMainVertex[til_]:=Map[add,til[[4,All,1]]]

GetPoints2[til_,F1_,F0_,theta_:2Pi/5]:=
	Module[{r=StarVectors[F1,F0,theta]},r[[1]]=Map[Normalize,r[[1]]];Map[add2[#,r]&,til[[2]],{2}]]
GetTiles2[til_,F1_,F0_,theta_:2Pi/5]:=
	Module[{r=StarVectors[F1,F0,theta]},r[[1]]=Map[Normalize,r[[1]]];Map[add2[#,r]&,til[[4]],{2}]]
GetLines2[til_,F1_,F0_,theta_:2Pi/5]:=
	Module[{r=StarVectors[F1,F0,theta]},r[[1]]=Map[Normalize,r[[1]]];Map[add2[#,r]&,til[[3]],{2}]]
GetMainVertex2[til_,F1_,F0_,theta_:2Pi/5]:=
	Module[{r=StarVectors[F1,F0,theta]},r[[1]]=Map[Normalize,r[[1]]];Map[add2[#,r]&,til[[4,All,1]]]]


GetPoints3[til_,r_]:=
	Module[{},Map[add2[#,r]&,til[[2]],{2}]]
GetTiles3[til_,angles_,theta_,indexlist_]:=
	Module[{r=StarVectors2[angles,theta,indexlist]},r[[1]]=Map[Normalize,r[[1]]];Map[add2[#,r]&,til[[4]],{2}]]
GetLines3[til_,r_]:=
	Module[{},Map[add2[#,r]&,til[[3]],{2}]]
GetMainVertex3[til_,angles_,theta_,indexlist_]:=
	Module[{r=StarVectors2[angles,theta,indexlist]},r[[1]]=Map[Normalize,r[[1]]];Map[add2[#,r]&,til[[4,All,1]]]]


(********************************************************************)
(*                      End of Private context                      *)
(********************************************************************)

End[]

(********************************************************************)
(*                                                                  *)
(*                    Protect exported functions                    *)
(*                                                                  *)
(********************************************************************)

Protect[StarVectors,
	PlotGridVectors,
	PlotGrid,
	GetIntersections,
	PlotIntersections,
	DualizeGrid,
	GetSingleTile,
	PlotDualTiling,
	PlotSingular,
	PlotColors,
	PlotColors2,
	GetPoints,
	GetPoints2,
	GetLines,
	GetTiles,
	GetMainVertex,
	GetLines2,
	GetTiles2,
	GetMainVertex2,
	GetIntersectionTiles,
	PlotIntersectionTiles,
	GetUnitCell,
	PlotUnitCell,
	add,
	add2,
	TrimUnitCell,
	GetAmmannUnitCell,
	FindNeighbors,
	PlotStarVectors,
	StarPhases2,
	StarVectors2,
	DualizeGrid2,
	PlotDualTiling2,
	PlotGrid2,
	PlotIntersections2,
	GetUnitCell2,
	TrimUnitCell2,
	TrimUnitCell3,
	PlotUnitCell2,
	PlotUnitCell3,
	GetUnitCell3,
	GetIndexList,
	GetKLattice,
	DualizeUnitCell,
	PlotDualUnitCell,
	PlotBZ,
	StarVectors3,
	GetInfoGraphics,
	ModStarVectors,
	GetEffectiveN,
	GetMinDist]


(********************************************************************)
(*                                                                  *)
(*        End of package "Tilings`GridMethod2DPA`"           *)
(*                                                                  *)
(********************************************************************)

EndPackage[]
