import "../BasicData.m": n;
import "../Initialize.m": K,tau,ZB,w,B,p1,p2;


matbas:=function(GG);
 //Internal method
 //Convert (2n)*(2n) matrice over Z into n*n matrices over O_K
 //Input&Output: List
 n:=NumberOfRows(GG[1]) div 2;
 MM:=[];
 for g in GG do
  M:=MatrixRing(K,n) ! 0;
  for i in [1..n-1] do
   for j in [1..n-1] do
    M[i][j]:= g[2*i-1][2*j-1] + g[2*i-1][2*j] * tau;
   end for;
   M[i][n]:= ZB[1]*g[2*i-1][2*n-1] + ZB[2]*g[2*i-1][2*n] ;
  end for;
  for j in [1..n-1] do
   M[n][j]:= g[2*n-1][2*j-1]/ZB[1]+tau*g[2*n-1][2*j]/ZB[1];
  end for;
  M[n][n]:=g[2*n-1][2*n-1]+ZB[2]/ZB[1]*g[2*n-1][2*n];
  Append(~MM,M);
 end for;
 return MM;
end function;

RealPart:=function(x)
 return (1/2)*(x+Conjugate(x));
end function;

ImaginaryPart:=function(x)
 return (1/(2*w))*(x-Conjugate(x));
end function;

matbas2:=function(L)
 //Internal method
 //Convert n*n matrices into (2n)*(2n) matrices over Z
 //[For matrix groups]
 res:=[];
 //Define Basis matrices for p1 and p2
 BM1:=MatrixRing(K,2)![[1,0],[RealPart(tau),ImaginaryPart(tau)]];
 BM2:=MatrixRing(K,2)![[RealPart(ZB[1]),ImaginaryPart(ZB[1])] ,
      [RealPart(ZB[2]) , ImaginaryPart(ZB[2])]];
 for x in L do
  M:=KMatrixSpace(K,0,2*n)![];
  //Compute images of basis vectors under x:
  ims:=[b*x : b in B];
  for i in [1..2*n] do
   //Compute their coefficients in the Z-basis:
   v:=ims[i][1];
   coeffs:=KMatrixSpace(K,1,0)![];
   for k in [1..n-1] do
    coeffs:=HorizontalJoin(coeffs,Solution(BM1,KMatrixSpace(K,1,2)![RealPart(v[k]),ImaginaryPart(v[k])]));
   end for;
   coeffs:=HorizontalJoin(coeffs,Solution(BM2,KMatrixSpace(K,1,2)![RealPart(v[n]),ImaginaryPart(v[n])]));
   M:=VerticalJoin(M,coeffs);
  end for;
  Append(~res,MatrixRing(Integers(),2*n)!M);
 end for;
 return res;
end function;

matbas3:=function(L)
 //Internal method
 //Convert n*n matrices into (2n)*(2n) matrices over Z
 //[For arbitrary matrices]
 res:=[];
 //Define Basis matrices for p1 and p2
 BM1:=MatrixRing(K,2)![[1,0],[RealPart(tau),ImaginaryPart(tau)]];
 //Use the standard vector space basis instead of the Z-basis for p2
 BM2:=BM1;
 for x in L do
  M:=KMatrixSpace(K,0,2*n)![];
  //Compute images of basis vectors under x:
  ims:=[b*x : b in B];
  for i in [1..2*n] do
   //Compute their coefficients in the Z-basis:
   v:=ims[i][1];
   coeffs:=KMatrixSpace(K,1,0)![];
   for k in [1..n-1] do
    coeffs:=HorizontalJoin(coeffs,Solution(BM1,KMatrixSpace(K,1,2)![RealPart(v[k]),ImaginaryPart(v[k])]));
   end for;
   coeffs:=HorizontalJoin(coeffs,Solution(BM2,KMatrixSpace(K,1,2)![RealPart(v[n]),ImaginaryPart(v[n])]));
   M:=VerticalJoin(M,coeffs);
  end for;
  //Here: MatrixRing(Rationals(),...) [for arbitrary matrices over K]
  Append(~res,MatrixRing(Rationals(),2*n)!M);
 end for;
 return res;
end function;

intrinsic ElementsOfNorm(norm::Any,P::RngQuadIdl) ->SeqEnum
 {A procedure to produce all elements of given norm in an ideal P over the integers of a number field}
 OK:=Integers(K);
 Gram:=MatrixRing(Integers(),2)![[2,Trace(tau)],[Trace(tau),2*Norm(tau)]];//We will use the trace form and use the magma procedures for Z-lattices
 L:=LatticeWithGram(Gram);
 S:=ShortVectors(L,2*norm,2*norm);//This is not exactly the ideal way since we first compute all vectors in O_K and then intersect with P
 output:=[s[1][1]+s[1][2]*tau: s in S];
 output:=[x: x in output | x in P];
 return output;
end intrinsic;


intrinsic ChooseBasis(L::SeqEnum,K::Fld)->SeqEnum
 {Choose a Basis of a RowSpace out of a given List of RowVectors via simple linear algebra}
 Indices:=[];
 T:=KMatrixSpace(K,#L,Dimension(Parent(L[1])))!L;
 T:=Transpose(T);
 T:=EchelonForm(T);
 for i in [1..NumberOfRows(T)] do
  for j in [1..NumberOfColumns(T)] do
   if T[i][j] ne 0 then	
    Append(~Indices,j); break;
   end if;
  end for;
 end for;
 return [L[i]: i in Indices];
end intrinsic;


intrinsic AllMinVecs(F::HermForm) -> SeqEnum
 {Take a list of vectors and return all multiples of them which still lie in the lattice and have the same idealnorm}
 if assigned F`AllMinVecs then
  return F`AllMinVecs;
 else
  S:=MinimalVectors((F));
  output:={s: s in S};
  if n eq 2 then 
   for v in S do 
    v10:=true;
    if v[1][1] ne 0 then
     X:=[w/v[1][1]: w in ElementsOfNorm(Norm(v[1][1]),p1)];//This is a list of elements of norm 1 which might fulfill xv \in L.
     v10:=false;
    end if;
    if v10 then
     X:=[w/v[1][2]: w in ElementsOfNorm(Norm(v[1][2]),p2)];
    else
     X:=[x:x in X| x*v[1][2] in p2];
    end if;
    output:=output join {x*v: x in X};
   end for;
   output:=output join {-x : x in output}; // This is probably not necessary 
  end if;
  if n eq 3 then
   for v in S do 
    v10:=true; 
    v20:=true;
    if v[1][1] ne 0 then
     X:=[w/v[1][1]: w in ElementsOfNorm(Norm(v[1][1]),p1)];//This is a list of elements of norm 1 which might fulfill xv \in L.
     v10:=false;
     v20:=false;
    end if;
    if v10 then
     if v[1][2] ne 0 then
      X:=[w/v[1][2]: w in ElementsOfNorm(Norm(v[1][2]),p1)];
      v20:=false;
     end if;
    end if;
    if v20 then 
     X:=[w/v[1][3]: w in ElementsOfNorm(Norm(v[1][3]),p2)];
    end if;
    X:=[x:x in X| x*v[1][3] in p2 and x*v[1][2] in p1 and x*v[1][1] in p1];
    output:=output join {x*v: x in X};
   end for;
   output:=output join {-x : x in output}; // This is probably not necessary 
  end if;
  F`AllMinVecs:=[v: v in output];
  return AllMinVecs(F);
 end if;
end intrinsic;


intrinsic LowIndexStabilizer(F::HermForm,CheckMembership::Any) -> Grp
 {A function to determine the stabilizer in a subgroup of GL which comes with a function CheckMembership (which does just that)}
 S:=StabilizerOfMinimalClass(F);
 Elements:=[x: x in S | CheckMembership(x)];
 order:=Order(MatrixGroup<n,K|Elements>);
 Gens:=[Elements[1]];
 H:=MatrixGroup<n,K|Elements>;
 while Order(H) ne order do
  x:=Random(Elements);
  if not x in H then
   Append(~Gens,x);
   H:=MatrixGroup<n,K|Gens>;
  end if;
 end while;
 return H;
end intrinsic;

U,u:=UnitGroup(Integers(K));
UU:={u(x): x in U};
trivialcenter:={MatrixRing(Integers(K),n)!u(x): x in U};

intrinsic BoundaryEmbeddings(LargeFace::HermForm,SmallFace::HermForm,CheckMembership::Any)->SeqEnum
 {Find Matrices in GL(L) such that gSmallFace is in the boundary of LargeFace}
 //Both LargeFace and SmallFace need to be given by a representative of the minimal class, CheckMembership will check wether a matrix is in GL(L), We want  to construct g in GL(L) such that minvecs(LargeFace)g subset minvecs(SmallFace)
 Target:=AllMinVecs(SmallFace);  //These are all possible images of minvecs(LargeFace) 
 M:=MinimalVectors((LargeFace));
 Source:=ChooseBasis(M,K); 
 Inv:=(MatrixAlgebra(K,n)!Source)^(-1);
 Invdet:=Determinant(Inv);
 ImageSets:=Subsets({x: x in Target},n);
 ImageSets:=[IndexedSetToSequence(SetToIndexedSet(x)): x in ImageSets];
 ImageSets:=[y: y in ImageSets | Invdet*Determinant(Matrix(y)) in UU];
 Sn:=SymmetricGroup(n);
 ImageLists:=[];
 for g in Sn do 
  for x in ImageSets do	
   Append(~ImageLists,[x[i^g]: i in [1..n]]);
  end for;
 end for;
 //Norms:=[IdealNorm(v): v in Source];
 //ImageLists:=[x:x in ImageLists|  Determinant(MatrixAlgebra(K,n)!x) in UU]; // Now ImageList is a list of all possible tuples [w1,...,wn] such that v1g=w1,...,vng=wn possibly determines an element g in GL(L) which might fit the bill
 //ImageLists:=[x: x in ImageLists| [IdealNorm(y): y in x] eq Norms]; // Since idealnorm(x)=idealnorm(xg) forall x and g in GL we need only consider lists which sequences of norms coincide with those of Source
 PossibleElementsList:=[Inv*(MatrixAlgebra(K,n)!x): x in ImageLists]; // This is a list of all possible GroupElements such that gSmallFace might be in the boundary of LargeFace
 PossibleElementsList:=[g: g in PossibleElementsList | CheckMembership(g)]; //Now only elements of GL(L) remain
 PossibleElementsList:=[g: g in PossibleElementsList | {m*g: m in M} subset {v: v in Target}]; // Now only elements remain which really fulfill gSmallFace in the boundary of LargeFace
 if #PossibleElementsList eq 0 then 
  return [];
 end if;
 FinalOutput:=[PossibleElementsList[1]];
//Now we will make the list duplicate free
 for g in PossibleElementsList do
  bool:=true;
  for h in FinalOutput do
   if h^(-1)*g in StabilizerOfMinimalClass(SmallFace) then	
    bool:=false;
    break h;
   end if;
  end for;
  if bool then	
   Append(~FinalOutput,g);
  end if;
 end for;
 return FinalOutput;
end intrinsic;



intrinsic Vertices(F::Mtrx,V::VData) -> SeqEnum
 {This will give a list of all perfect forms in the boundary of a minimal class represented by F}
 if not assigned(V`MinimalClasses) then
  MinimalClasses(V);
 end if;
 FormREPS:=V`MinimalClasses;
 output:=[**];
 for i in [1..#FormREPS[1]] do
  for x in BoundaryEmbeddings(F,FormREPS[1][i]`Matrix, IsInGL) do
   Append(~output,[*i,x*]);
  end for;
 end for;
 return output;
end intrinsic;


U,u:=UnitGroup(Integers(K));
trivialcenter:={MatrixRing(Integers(K),n)!u(x): x in U};


intrinsic OrientationSignByDeterminant(F::HermForm,g::Mtrx,W::VData) -> RngIntElt 
 {Computes the orientation action of g on the minimal class represented by F via the use of determinants}
 if g in trivialcenter then
  return 1;
 end if;
 if not assigned(W`MinimalClasses) then
  MinimalClasses(W);
 end if;
 FormREPS:=W`MinimalClasses;
 
 
 //Let's try something else (instead of vertices)
 k:=PerfectionCorank(F);
 if k eq 0 then
  return 1;
 end if;
 if assigned F`TSpaceBasis and assigned F`TSpaceRat then
  TSpaceBasis:=F`TSpaceBasis;
  Space:=F`TSpaceRat;
 else
  output:=[**];
  kk:=k eq 0 select 0 else k-1;
  for i in [1..#FormREPS[kk+1]] do
   for x in BoundaryEmbeddings(F,FormREPS[kk+1][i], IsInGL) do
    Append(~output,[*i,x*]);
   end for;
  end for;
  V:=output;
  V:=[y[2]*FormREPS[kk+1][y[1]]`Matrix*HermitianTranspose(y[2]): y in V];
  //Until here is something else
  //V:=[y[2]*FormREPS[1][y[1]]`Matrix*HermitianTranspose(y[2]): y in V];
  GensOfTSpace:=[V[1]-V[i]: i in [1..#V]];
  GensOfTSpaceRat:=[ElementToSequence(y): y in matbas3(GensOfTSpace)];
  Indices:=[];         //Now we choose a (Q)-Basis for the Space of Translations of affine space generated by the given class
  T:=KMatrixSpace(Rationals(),#GensOfTSpaceRat,4*n^2)!GensOfTSpaceRat;
  assert Rank(T) eq k;
  T:=Transpose(T);
  T:=EchelonForm(T);
  for i in [1..NumberOfRows(T)] do
   for j in [1..NumberOfColumns(T)] do
    if T[i][j] ne 0 then
     Append(~Indices,j); break;
    end if;
   end for;
  end for;
  F`TSpaceBasis:=[GensOfTSpace[i]: i in Indices];
  TSpaceBasis:=F`TSpaceBasis;
  F`TSpaceRat:=KSpaceWithBasis(KMatrixSpace(Rationals(),#TSpaceBasis,4*n^2)![ElementToSequence(y): y in matbas3(TSpaceBasis)]);
  Space:=F`TSpaceRat;
 end if;
 Images:=[g*TSpaceBasis[i]*HermitianTranspose(g): i in [1..#TSpaceBasis]];
 Images:=matbas3(Images);
 Images:=[ElementToSequence(y): y in Images];
 //Space:=KSpaceWithBasis(KMatrixSpace(Rationals(),#TSpaceBasis,4*n^2)![ElementToSequence(y): y in matbas3(TSpaceBasis)]);
 MatrixRep:=[Coordinates(Space,Space!y):y in Images];   //We determine the basis representation of the images under g
 return Sign(Determinant(Matrix(MatrixRep)));
end intrinsic;

intrinsic EvenStabilizerOfMinimalClass(F::HermForm,V::VData) -> Grp
 {The orientation preserving subgroup of the stabilizer of the minimal class}
 if assigned F`EvenMinClassStabilizer then
  return F`EvenMinClassStabilizer;
 end if;
 Stab:=StabilizerOfMinimalClass(F);
 if PerfectionCorank(F) eq 0 or #Stab eq #trivialcenter then
  return Stab;
 end if;
 StabGen:=[x: x in Generators(Stab)];
 Img:=[];
 for x in StabGen do
  Append(~Img,<x,GL(1,Integers())![OrientationSignByDeterminant(F,x,V)]>);
 end for;
 orientationchar:=hom<Stab->GL(1,Integers())|Img>;
 F`EvenMinClassStabilizer:=Kernel(orientationchar);
 return F`EvenMinClassStabilizer;
end intrinsic;


 
intrinsic EvenStabilizerOfMinimalClass2(F::HermForm,V::VData) -> Grp
 {This will compute the orientation preserving subgroup of the stabilizer of a minimal class, this is a normal subgroup of index 1 or 2 }
 Stab:=StabilizerOfMinimalClass(F);
 index:=1;
 EvenGenerators:={};
 for i in [1..#Generators(Stab)] do  //Let us first determine whether the index is 1 or 2.
   g:=Stab.i;             
  if OrientationSignByDeterminant(F,g,V) eq -1 then
   index:=2;
   Include(~EvenGenerators,g^2);
   ii:=i;
   break;
  else
   Include(~EvenGenerators,g);
  end if;
 end for;
 if index eq 1 then 
  return Stab;
 end if;
 Ordnung:=Order(Stab) div 2;
 S:=MatrixGroup<n,K|[y : y in EvenGenerators]>;
 if Order(S) eq Ordnung then
  return S;
 end if;
 for i in [ii+1..#Generators(Stab)] do
  g:=Stab.i;
  if OrientationSignByDeterminant(F,g,V) eq -1 then
   Include(~EvenGenerators,g^2);
  else
   Include(~EvenGenerators,g);
  end if;
  S:=MatrixGroup<n,K|[y : y in EvenGenerators]>;
  if Order(S) eq Ordnung then
   return S;
  end if;
 end for;

 while Order(S) ne Ordnung do    //We add new elements until the order is large enough
  x:=Random(Stab);
  if x in S then continue; end if;
  if OrientationSignByDeterminant(F,x,V) eq 1 then
   Include(~EvenGenerators,x);
  else
   Include(~EvenGenerators,x^2);
  end if;
  S:=MatrixGroup<n,K|[y : y in EvenGenerators]>;
 end while;
 return S;
end intrinsic;


intrinsic LowIndexEvenStabilizer(F::HermForm,CheckMembership::Any,V::VData) -> SeqEnum
 {This will compute the orientation preserving subgroup of the stabilizer of a minimal class (in a subgroup of GL), this is a normal subgroup of index 1 or 2 }
 Stab:=LowIndexStabilizer(F,CheckMembership);
 FullEvenStab:=EvenStabilizerOfMinimalClass(F,V);
 return Stab meet FullEvenStab;
end intrinsic;


intrinsic IsInSL(g::Mtrx) -> BoolElt
 {A function to check wether something is in SL}
 return IsInGL(g) and Determinant(g) eq 1;
end intrinsic;

//An implementation of product replacement:
Step:=function(T)
 L:=#(T);
 i:=Random(1,L);
 j:=Random(1,L);
 while i eq j 
  do j:=Random(1,L); 
 end while;
 eps:=Random(0,1);
 if eps eq 1 
  then T[i]:=T[i]*T[j]; 
 else 
  T[i]:=T[j]^(-1)*T[i]; 
 end if;
 return T;
end function;


PR:=function(X,B)
 T:=[];
 A:=X[1]*X[1]^(-1);
 for i in [1..Max(#(X)+2,10)] do
  if i le #(X) 
   then T[i]:=X[i]; 
  else 
   T[i]:=1; 
  end if; 
 end for;
 for i in [1..B] do 
  j:=Random(1,#(T)); 
  T:=Step(T); 
  eps:=Random(0,1);
  if eps eq 1 
   then A:=T[j]*A; 
  else 
   A:=T[j]^(-1)*A; 
  end if; 
 end for;
 return A;
end function;

 
intrinsic IsInCongruenceSubgroupTest(p::RngQuadIdl) ->Any
 {A function to produce the checkmembershipfunction for the full congruence subgroup of level p}
 f:=function(M)
  return IsInSL(M) and {M[1][1]-1 in p,M[1][2] in p, M[2][1] in p, M[2][2]-1 in p} eq {true};
 end function;
 return f;
end intrinsic;


intrinsic SystemOfRepresentativesFiniteIndex(gens::Any,CheckMembership::Any,ind::RngIntElt)-> SeqEnum
 {A function to compute a system of representatives in GL if we know the index of the subgroup and a function to check for membership}
 reps:=[Parent(gens[1])!1];
 iterations:=2;
 counter:=1;
 while #reps ne ind do
  if counter mod 100 eq 0 then
   iterations+:=1;
  end if;
  x:=PR(gens,iterations);
  if {CheckMembership(y*x^(-1)):y in reps} eq {false} then
   Append(~reps,x);
  end if;
 counter:=counter+1;
 end while;
 return reps;
end intrinsic;













