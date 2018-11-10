import "../../BasicData.m": n,d,steinitz;
import "../../Initialize.m": K,w,tau,p1,p2,ZB,B,BasHermNorm,mmax,F,sqrtd,Injec;
import "../Functions/Groups.m": matbas,matbas2,matbas3;

DimSym:=n^2;
ERank:=DimSym-2;
FRank:=DimSym-1;

function MatrixToCoordinates(A)
 coords:=[];
 for i in [1..n] do
  Append(~coords,Rationals()!A[i][i]);
  for j in [i+1..n] do
   coords cat:= Eltseq(A[i][j]);
  end for;
 end for;
 return coords;
end function;
 

intrinsic ComputeEdges(V::VData)
 {Compute and save the attribute Edges in V.}
 Edges:=[* *];
 FL:=V`FacesList;
 for i in [1..#FL] do
  edges:={* *};
  F:=FL[i];
  S:=(V`PerfectList[i])`Vertices; 
  for j in [1..#F] do
   for k in [j+1..#F] do
    s:=F[j] meet F[k];
    if s ne {} then
     s:=Setseq(F[j] meet F[k]);
     d:=PerfectionRankListOfMatrices([S[l] : l in s]);
     if d eq ERank then
      Include(~edges,s);
     end if;
    end if;
   end for;
  end for;
  Append(~Edges,edges);
 end for;
 V`Edges:=Edges;
end intrinsic;

intrinsic ComputeBarycenters(V::VData)
 {Compute and save the attribute Barycenters in V.}
 Barycenters:=[* *];
 PL:=V`PerfectList;
 for P in PL do
  S:=MinimalVectors(P);
  B:=(1/#S)*(&+[HermitianTranspose(x)*x/Trace(HermitianTranspose(x)*x) : x in S]);
  Append(~Barycenters,B);
 end for;
 V`Barycenters:=Barycenters;
end intrinsic;

intrinsic ComputeStabilizers(V::VData: SL:=false,CheckMembership:=0)
 {Compute and save the attribute Stabilizers in V.}
 PL:=V`PerfectList;
 V`Stabilizers:=[* ConvertGroupToNumberField(HermitianAutomorphismGroup(P:SL:=SL,CheckMembership:=CheckMembership)) : P in PL *];
end intrinsic;

intrinsic GeodesicIntersection(M1::Mtrx,M2::Mtrx,Face::Any : printsolution:=false) -> BoolElt
 {Checks if the geodesic from M1 to M2 intersects the face "Face".}
 Q:=Rationals();
 M1:=M1/Trace(M1);
 M2:=M2/Trace(M2);
 Face:=[x/Trace(x) : x in Face];
 Translations:=[Face[i]-Face[1] : i in [2..#Face]];
 A:=KMatrixSpace(Rationals(),#Face,n^2) ! 0;
 for i in [1..#Translations] do
  a:=[];
  for k in [1..n] do
   Append(~a,Rationals()!Translations[i][k][k]);
   for l in [k+1..n] do
    a cat:= Eltseq(Translations[i][k][l]);
   end for;
  end for;
  for j in [1..#a] do
   A[i][j]:=a[j];
  end for;
 end for;
 a:=[];
 for i in [1..n] do
  Append(~a,Rationals()!( (M1-M2)[i][i] ) );
  for j in [i+1..n] do
   a cat:=Eltseq((M1-M2)[i][j]);
  end for;
 end for;
 for i in [1..#a] do
  A[#Face][i]:=a[i];
 end for;
 LS:=[];
 for i in [1..n] do
  Append(~LS,Rationals()!( (M1-Face[1])[i][i] ) );
  for j in [i+1..n] do
   LS cat:=Eltseq((M1-Face[1])[i][j]);
  end for;
 end for;
 LS:=Vector(LS);
 a,b,c:=IsConsistent(A,LS);
 if not a then 
  return false,false;
 end if;
 t:=b[#Face];
 if not (0 lt t and t lt 1) then
  return false,false;
 end if;
 GeneratorsOfPolytope:=[];
 for f in Face do
  coords:=[];
  for i in [1..n] do
   Append(~coords,Rationals()!f[i][i]);
   for j in [i+1..n] do
    coords cat:= Eltseq(f[i][j]);
   end for;
  end for;
  Append(~GeneratorsOfPolytope,coords);
 end for;
 P:=Polytope(GeneratorsOfPolytope);
 solution:=M1+t*(M2-M1);
 coords:=[];
 for i in [1..n] do
  Append(~coords,Q!solution[i][i]);
  for j in [i+1..n] do
   coords cat:= Eltseq(solution[i][j]);
  end for;
 end for;
 Q:=Polytope([coords]);
 if printsolution then
  result:=MatrixRing(K,n)!0;
  k:=0;
  for i in [1..n] do
   result[i][i]:=coords[1+k]; k+:=1;
   for j in [i+1..n] do
    result[i][j]:=coords[k+1]+w*coords[k+2]; k+:=2;
    result[j][i]:=ComplexConjugate(result[i][j]);
   end for;
  end for;
  return Q subset P, result;
 end if;
 return Q subset P , false;
end intrinsic;

intrinsic GeodesicIntersection2(M1::Mtrx,M2::Mtrx,Face::Any : printsolution:=false) -> BoolElt
 {Checks if the geodesic from M1 to M2 intersects the face "Face".}
 Q:=Rationals();
 M1:=M1/Trace(M1);
 M2:=M2/Trace(M2);
 Face:=[x/Trace(x) : x in Face];
 Translations:=[Face[i]-Face[1] : i in [2..#Face]];
 A:=KMatrixSpace(Rationals(),#Face,n^2) ! 0;
 for i in [1..#Translations] do
  a:=[];
  for k in [1..n] do
   Append(~a,Rationals()!Translations[i][k][k]);
   for l in [k+1..n] do
    a cat:= Eltseq(Translations[i][k][l]);
   end for;
  end for;
  for j in [1..#a] do
   A[i][j]:=a[j];
  end for;
 end for;
 a:=[];
 for i in [1..n] do
  Append(~a,Rationals()!( (M1-M2)[i][i] ) );
  for j in [i+1..n] do
   a cat:=Eltseq((M1-M2)[i][j]);
  end for;
 end for;
 for i in [1..#a] do
  A[#Face][i]:=a[i];
 end for;
 LS:=[];
 for i in [1..n] do
  Append(~LS,Rationals()!( (M1-Face[1])[i][i] ) );
  for j in [i+1..n] do
   LS cat:=Eltseq((M1-Face[1])[i][j]);
  end for;
 end for;
 LS:=Vector(LS);
 a,b,c:=IsConsistent(A,LS);
 if not a then 
  return false,false;
 end if;
 t:=b[#Face];
 if not (0 le t and t le 1) then
  return false,false;
 end if;
 GeneratorsOfPolytope:=[];
 for f in Face do
  coords:=[];
  for i in [1..n] do
   Append(~coords,Rationals()!f[i][i]);
   for j in [i+1..n] do
    coords cat:= Eltseq(f[i][j]);
   end for;
  end for;
  Append(~GeneratorsOfPolytope,coords);
 end for;
 P:=Polytope(GeneratorsOfPolytope);
 solution:=M1+t*(M2-M1);
 coords:=[];
 for i in [1..n] do
  Append(~coords,Q!solution[i][i]);
  for j in [i+1..n] do
   coords cat:= Eltseq(solution[i][j]);
  end for;
 end for;
 Q:=Polytope([coords]);
 if printsolution then
  result:=MatrixRing(K,n)!0;
  k:=0;
  for i in [1..n] do
   result[i][i]:=coords[1+k]; k+:=1;
   for j in [i+1..n] do
    result[i][j]:=coords[k+1]+w*coords[k+2]; k+:=2;
    result[j][i]:=ComplexConjugate(result[i][j]);
   end for;
  end for;
  return Q subset P, result;
 end if;
 return Q subset P , false;
end intrinsic;

intrinsic ConstructWordViaBarycenters(x::Mtrx,V::VData) -> SeqEnum
 {Solve the word problem for x in GL(L) using the Voronoi data V using the barycenter method. This will not necessarily terminate.}
 if not assigned V`Stabilizers then ComputeStabilizers(V); end if;
 if not assigned V`Barycenters then ComputeBarycenters(V); end if;
 Barycenters:=V`Barycenters;
 Stabs:=&cat[  [g : g in G] : G in V`Stabilizers];
 PL:=V`PerfectList;
 FL:=V`FacesList;
 FTL:=V`FaceTrafoList;
 NL:=V`NeighbourList;
 OKGens:=V`OKGens;
 B1:=Barycenters[1];
 B:=B1;
 B2:=HermitianTranspose(x)*B1*x;
 Currents:=[1];
 cur:=1; //current index
 Word:=[];
 found:=true;
 while (not (x in Stabs)) and found do
  found:=false;
  for j in [1..#FL[cur]] do
   S:=MinimalVectors(PL[cur]);
   FF:=[HermitianTranspose(S[i])*S[i] : i in FL[cur][j]];
   if GeodesicIntersection(B1,B2,FF) then
    found:=true;
    g:=FTL[cur][j];
    cur:=NL[cur][j]; Append(~Currents,cur);
    x:=x*g^-1; 
    if g ne MatrixRing(K,n)!1 then  
     p:=0; e:=0;
     if g in OKGens then p:=Position(OKGens,g); e:=1; end if;
     if -g in OKGens then p:=Position(OKGens,-g); e:=1; end if;
     if g^-1 in OKGens then p:=Position(OKGens,g^-1); e:=-1; end if;
     if -g^-1 in OKGens then p:=Position(OKGens,-g^-1); e:=-1; end if;
     if [p,e] eq [0,0] then return cur,j; end if;
     Append(~Word,[p,e]);
    end if;
    B1:=Barycenters[cur];
    B2:=HermitianTranspose(x)*B*x;
    break;
   end if;
  end for;
 end while;
 Append(~Word,[0,Position(Stabs,x)]);
 return found,Word,x,Currents,B;
end intrinsic;

intrinsic ConstructWord(x::Mtrx,V::VData:SL:=false,CheckMembership:=0) -> SeqEnum
 {Solve the word problem for x in GL(L) by the method described in the paper using the Voronoi data V }
 if not assigned V`Stabilizers then ComputeStabilizers(V:SL:=SL,CheckMembership:=CheckMembership); end if;
 if not assigned V`Barycenters then ComputeBarycenters(V); end if;
 if not assigned V`Edges then ComputeEdges(V); end if;
 Edges:=V`Edges;
 Barycenters:=V`Barycenters;
 Stabs:=[g : g in V`Stabilizers[1]];
 PL:=V`PerfectList;
 FL:=V`FacesList;
 FTL:=V`FaceTrafoList;
 NL:=V`NeighbourList;
 OKGens:=V`MultFreeList;
 xoriginal:=x;
 B1:=Barycenters[1];
 B:=B1;
 B2:=HermitianTranspose(x)*B1*x;
 cur:=1; //current index
 S:=MinimalVectors(PL[1]);
 Word:=[**];
 found:=true;
 while not x in V`Stabilizers[1] do
  if not found then
   for indices in Edges[cur] do
    if GeodesicIntersection2(B1,B2,[HermitianTranspose(S[i])*S[i] : i in indices]) then
     B1:=Barycenters[cur];
     changed:=false;
     M:=MinimalVectors(PL[cur]);
     MM:=[HermitianTranspose(x)*x/Trace(HermitianTranspose(x)*x) : x in M];
     N:=1;
     while not changed do
      N*:=100;
      B1pert:=B1;
      for i in [1..n] do
       B1pert[i][i]+:=1/Random(N,2*N);
       for j in [i+1..n] do
        r1:=1/Random(N,2*N); r2:=1/Random(N,2*N);
        B1pert[i][j]+:=r1+w*r2; B1pert[j][i]+:=r1-w*r2;
       end for;
      end for;
      B1pert/:=Trace(B1pert);
      Q:=Polytope([ MatrixToCoordinates(B1pert) ]);
      P:=Polytope([ MatrixToCoordinates(x) : x in MM ]);
      if Q subset P then
       changed:=true;
       B1:=B1pert;
      end if;
     end while;
     break indices;
    end if;
   end for;
  end if;
  found:=false;
  S:=MinimalVectors(PL[cur]);
  for j in [1..#FL[cur]] do
   FF:=[HermitianTranspose(S[i])*S[i] : i in FL[cur][j]];
   b,r:=GeodesicIntersection(B1,B2,FF:printsolution:=true);
   if b then
    found:=true; 
    g:=FTL[cur][j];
    cur:=NL[cur][j];
    x:=x*g^-1; 
    if g ne MatrixRing(K,n)!1 then  
     p:=0; e:=0;
     if g in OKGens then p:=Position(OKGens,g); e:=-1; end if;
     if [p,e] eq [0,0] then return cur,j; end if;
     Append(~Word,[p,e]);
    end if;
    B1:=HermitianTranspose(g^-1)*r*g^-1;
    B2:=HermitianTranspose(x)*B*x;
    break;
   end if;
  end for;
 end while;
 Append(~Word,[0,Position(Stabs,x)]);
 //At this point the output is to be read as follows:
 //x*w[1]^exp[1]*w[2]^exp[2]*....*w[r]^exp[r] = Stabs[l]
 //The stabilizer element is found in the last entry of the list "Word".

 //For convenience, we change this, so that now we have x=s*w[2]^exp[2]*...*w[r]^exp[r] and s is the stabilizer element, which is now found in the first entry
 temp:=[ Word[#Word] ];
 for i in [1..#Word-1] do
  Word[#Word-i][2]*:=-1;
  Append(~temp,Word[#Word-i]);
 end for;
 Word:=temp;
 return found,Word,x;
end intrinsic;


intrinsic SolveWordProblem(x::Mtrx,V::VData)->GrpFpElt
 {}
 found,Word,xx:=ConstructWord(x,V);
 if not found then
  error "This element does not belong to the computed group";
 end if;
 Stabs:=[g : g in V`Stabilizers[1]];

 G,mathom:=Presentation(V);
 simpl:=V`GNonsimplifiedToSimplifiedHom;
 Gnon:=V`GNonsimplified;
 res:=G!1;
 for i in [2..#Word] do
  res:=res*simpl(Gnon.Word[i][1]^Word[i][2]);
 end for;

 s:=Stabs[Word[1][2]];
 pos:=Position([P[1]: P in V`StabilizerPreimages],s);
 res:=simpl(V`StabilizerPreimages[pos][2])*res;
 assert V`GSimplifiedHom(res) eq x;
 
 return res;
end intrinsic;