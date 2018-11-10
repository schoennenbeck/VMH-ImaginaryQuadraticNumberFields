import "../../BasicData.m": n,d,steinitz;
import "../../Initialize.m": K,w,tau,p1,p2,ZB,B,BasHermNorm,mmax,F,sqrtd,Injec;
import "../Functions/Groups.m": matbas,matbas2,matbas3;

intrinsic ComputeRelations(V::VData:SL:=false,CheckMembership:=0)
 {Compute the relations corresponding to "circling" edges w.r.t. V.}
 if not assigned V`Stabilizers then ComputeStabilizers(V:SL:=SL,CheckMembership:=CheckMembership); end if;
 if not assigned V`Barycenters then ComputeBarycenters(V); end if;
 if not assigned V`Edges then ComputeEdges(V); end if;
 Barycenters:=V`Barycenters;
 Edges:=V`Edges;
 Stabs:=&cat[  [g : g in G] : G in V`Stabilizers];
 PL:=V`PerfectList;
 FL:=V`FacesList;
 FTL:=V`FaceTrafoList;
 NL:=V`NeighbourList;
 OKGens:=V`MultFreeList;
 Words:=[];
 Relations:=[];
 for a in [1..#PL] do
  for E in Edges[a] do
   for f in [1..#FL[a]] do
    if E subset FL[a][f] then
     cur:=a;
     S:=PL[cur]`Vertices;
     face:=f;
     groupelt:=FTL[cur][face];
     groupeltprevious:=Parent(groupelt)!1;
      word:=[[cur,face]];
      p:=0; e:=0;
      if groupelt in OKGens then p:=Position(OKGens,groupelt); e:=1; end if; 
      relation:=p eq 0 select [**] else [*<p,e>*];
      Nb:=NL[cur][face];
      finished:=false;
      k:=1;
      count:=1;
      while not finished do
       count +:= 1; 
       h1:=groupeltprevious^-1;
       F:={HermitianTranspose(h1^-1)*PL[cur]`Vertices[x]*h1^-1 : x in FL[cur][face]};  //!!!
       h:=groupelt^-1;
       X:=[{HermitianTranspose(h^-1)*PL[Nb]`Vertices[x]*h^-1 : x in f} : f in FL[Nb]];
       for i in [1..#X] do
        if {S[i] : i in E} subset X[i] and X[i] ne F then
         cur:=Nb;
         face:=i;
         g:=FTL[cur][face];
         if g ne MatrixRing(K,n)!1 then
          Append(~word,[Nb,i]);
          p:=0; e:=0;
          if g in OKGens then p:=Position(OKGens,g); e:=1; end if;          
          Append(~relation,<p,e>);
         end if;
         Nb:=NL[cur][face];
         groupeltprevious:=groupelt;
         groupelt:=g*groupelt;
         break;
        end if;
       end for;
       k+:=1;      
       if groupelt in V`Stabilizers[a] then
        Append(~relation,<a,groupelt>);
        finished:=true;
       end if;
      end while;
      Append(~Words,word);
      Append(~Relations,relation);
    end if;
   end for;
  end for;
 end for;
 V`Words:=Words;
 V`Relations:=Relations;
end intrinsic;




 
 


intrinsic Presentation(V::VData: projective:=false, simplify:=true,SL:=false,CheckMembership:=0)-> GrpFP
 {Computes a presentation of the arithmetic group underlying V}
 if assigned V`Presentation and V`Presentation then
  if projective then
   if simplify then
    return V`GSimplified,V`GSimplifiedHom,V`GSimplifiedPro,V`GSimplifiedToProHom;
   else
    return V`GNonsimplified, V`GNonsimplifiedHom,V`GNonsimplifiedPro,V`GNonsimplifiedToProHom;
   end if;
  else
   if simplify then
    return V`GSimplified,V`GSimplifiedHom;
   else
    return V`GNonsimplified, V`GNonsimplifiedHom;
   end if;
  end if;
 end if;
 if not assigned V`Relations then
  ComputeRelations(V:SL:=SL,CheckMembership:=CheckMembership);
 end if;
 Words:=V`Words;
 Gen:=V`MultFreeList;
 Rel:=V`Relations;
 Stabs:=V`Stabilizers;
 FPStabs:=[];
 FPStabsHom:=[**];
 PL:=V`PerfectList;
 FTL:=V`FaceTrafoList;
 FL:=V`FacesList;
 PNL:=V`PerfectNeighbourList;
 NL:=V`NeighbourList;
 One:=Parent(Gen[1])!1;
 for S in Stabs do
  F,f:=FPGroup(S);
  Append(~FPStabs,F);
  Append(~FPStabsHom,f);
 end for;
 ng:=#Gen+&+[#Generators(x): x in FPStabs];
 F:=FreeGroup(ng);
 R:=[];
 //print "A";
 //First compute the relations coming from the intersection of stabilizers
 for i in [1..#Stabs] do 
  for j in [i+1..#Stabs] do
   offseti:=i eq 1 select #Gen else #Gen+ &+[Ngens(FPStabs[k]): k in [1..i-1]];
   offsetj:=j eq 1 select #Gen else #Gen+ &+[Ngens(FPStabs[k]): k in [1..j-1]];
   for x in Stabs[i] meet Stabs[j] do
    if x ne Stabs[i]!1 then
     X:=x @@FPStabsHom[i];
     Y:=x @@ FPStabsHom[j];
     FX:=[k gt 0 select k+offseti else k-offseti : k in Eltseq(X)];
     FY:=[k gt 0 select k+offsetj else k-offsetj : k in Eltseq(Y)];
     Append(~R,F!FX=F!FY);
    end if;
   end for;
  end for;
 end for;
 
//print "B";

 //Now compute the relations due to the edges
 for w in Rel do
  lhs:=F!1;
  lhsmat:=One;
  for i in [1..#w-1] do
   lhs:=(F.w[i][1])^(w[i][2])*lhs;
   lhsmat:=(Gen[w[i][1]])^(w[i][2])*lhsmat;
  end for;
  cur:=w[#w][1];
  rhsmat:=w[#w][2];
  if rhsmat eq -lhsmat then
   rhsmat:=-rhsmat;
  end if;
  rhs:=rhsmat @@ FPStabsHom[cur];
  rhs:=ElementToSequence(rhs);
  offset:=cur eq 1 select #Gen else #Gen+&+[#Generators(FPStabs[i]): i in [1..cur-1]];
  rhs:=[Sign(i) gt 0 select i+offset else i-offset: i in rhs];
  Append(~R,lhs=F!rhs);
 end for;
 //print "C";
 //These are the relations coming from the finite stabilizers
 for cur in [1..#FPStabs] do 
  offset:=cur eq 1 select #Gen else #Gen+&+[#Generators(FPStabs[i]): i in [1..cur-1]];
  SRel:=Relations(FPStabs[cur]);
  for r in SRel do
   lhs:=ElementToSequence(LHS(r));
   rhs:=ElementToSequence(RHS(r));
   rhs:=[Sign(i) gt 0 select i+offset else i-offset: i in rhs];
   lhs:=[Sign(i) gt 0 select i+offset else i-offset: i in lhs];
   Append(~R,F!rhs=F!lhs);
  end for;
 end for;
//print "D";
 //These are the Poincare relations due to facets which lie in the same orbit.
 for p in [1..#PL] do
  for f in [1..#FL[p]] do
    nb:=NL[p][f];
    NB:=PL[nb];
    x:=FTL[p][f];
    X:=Position(Gen,FTL[p][f]);
    face:={PL[p]`Vertices[i]: i in FL[p][f]};
    Xinvface:={HermitianTranspose(x^-1)*v*x^-1: v in face};
    NBfaces:=[{NB`Vertices[i] : i in y}: y in FL[nb]];
    cor:=Position(NBfaces,Xinvface);
    y:=FTL[nb][cor];
    Y:=Position(Gen,y);
    lhs1:=[Y,X];    
    rhs1:=ElementToSequence((y*x) @@ FPStabsHom[p]);    
    offset1:=p eq 1 select #Gen else #Gen+&+[#Generators(FPStabs[i]): i in [1..p-1]];
    rhs1:=[Sign(i) gt 0 select i+offset1 else i-offset1: i in rhs1];  
    Append(~R,F!lhs1=F!rhs1);
   
   
  end for;
 end for;
//print "E";
 //Stabilizer Orbits of FTL
 for p in [1..#PL] do
  for f in [1..#FL[p]] do
   for f1 in [f+1..#FL[p]] do
    if NL[p][f] eq NL[p][f1] then
     x:=FTL[p][f];
     X:=Position(Gen,x);
     y:=FTL[p][f1];
     Y:=Position(Gen,y);
     for g in Stabs[p] do
      s:=x*g*y^-1;
      if s in Stabs[p] then
       offset:=p eq 1 select #Gen else #Gen+&+[#Generators(FPStabs[i]): i in [1..p-1]];
       lhs:=[X] cat [t gt 0 select t+offset else t-offset: t in Eltseq(g @@ FPStabsHom[p])] cat [-Y];
       rhs:=[t gt 0 select t+offset else t-offset: t in Eltseq(s @@ FPStabsHom[p])];
       Append(~R,F!lhs=F!rhs);
       break g;
      end if;
     end for;
    end if;
   end for;
  end for;
 end for;     
      
      
    //print "F";

 
 Img:=Gen;
 for i in [1..#PL] do
  Img cat:=[FPStabsHom[i](FPStabs[i].j): j in [1..#Generators(FPStabs[i])]];
 end for;

 GenLin:=sub<GeneralLinearGroup(n,K)|Img>;
 require IsSatisfied(R,[GeneralLinearGroup(n,K)!x: x in Img]): "Some of the computed relations are not satisfied. This should not have happened";
 G:=quo<F|R>;
 V`GNonsimplified:=G; 

 g:=hom<G->GL(n,K)|Img>;
 V`GNonsimplifiedHom:=g;
 
 H,h:=Simplify(G);
 hnew:=hom<H->GenLin|[g(h(H.i)): i in [1..#Generators(H)]]>;
 V`GSimplified:=H;
 V`GSimplifiedHom:=hnew; 
 V`GNonsimplifiedToSimplifiedHom:=h; 

 
 offset:=#Gen;
 if projective then
  NegOne:=(Parent(Gen[1])!-1) @@ FPStabsHom[1];
  NegOne:=Eltseq(NegOne);
  NegOne:=[t gt 0 select t + offset else t-offset: t in NegOne];
  RP:=R cat [F!NegOne=F!1];
  Gmod:=quo<F|RP>;
  pi:=hom<G->Gmod|[Gmod.i: i in [1..#Generators(Gmod)]]>;
  V`GNonsimplifiedPro:=Gmod;
  V`GNonsimplifiedToProHom:=pi;
  Hmod,hmod:=Simplify(Gmod);
  pinew:=hom<H->Hmod|[pi(h(H.i)) @@ hmod: i in [1..#Generators(H)]]>;
  V`GSimplifiedPro:=Hmod;
  V`GNonsimplifiedProToSimplifiedProHom:=hmod;
  V`GSimplifiedToProHom:=pinew;
 end if;


 

  //We now compute for each element in one of the finite stabilizers a preimage in the FPgrp
 StabsToG:=[**];
 for cur in [1..#FPStabs] do 
  offset:=cur eq 1 select #Gen else #Gen+&+[#Generators(FPStabs[i]): i in [1..cur-1]];
  for x in Stabs[cur] do
   w:=x @@FPStabsHom[cur];
   w:=Eltseq(w);
   w:=[t gt 0 select t + offset else t-offset: t in w];
   //g(G!w) eq x;
   Append(~StabsToG,<x,G!w>);
  end for;
 end for;

 V`StabilizerPreimages:=StabsToG;
 V`Presentation:=true; 
 return Presentation(V:projective:=projective,simplify:=simplify); 
end intrinsic;
 
