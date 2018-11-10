import "../../BasicData.m": n,d,steinitz;
import "../../Initialize.m": K,w,tau,p1,p2,ZB,B,BasHermNorm,mmax,F,sqrtd;

intrinsic HermitianAutomorphismGroup(A::HermForm : SL:=false, CheckMembership:=0) -> Grp
 {Returns the automorphism group of the input matrix A on the lattice O_K^(n-1)+p2 pecified in the BasicData of the package.}
 if CheckMembership cmpeq 0 then 
  if assigned A`Stabilizer and not SL then
   return A`Stabilizer;
  elif assigned A`SLStabilizer and SL then
   return A`SLStabilizer;
  elif not SL then   
   Ae:=TraceForm(A);
   mul:=Lcm([Denominator(x) : x in Eltseq(Ae)]);
   Ae:=ChangeRing(mul*Ae,Integers());
 
   Ae2:=TraceForm(NewHermForm(tau*A`Matrix));
   mul:=Lcm([Denominator(x) : x in Eltseq(Ae2)]);
   Ae2:=ChangeRing(mul*Ae2,Integers());
 
   L:=LatticeWithGram(Ae);
 
   G:=AutomorphismGroup(L,[Ae2]);
   A`Stabilizer:=G;
 
   return HermitianAutomorphismGroup(A);
  else
   G:=HermitianAutomorphismGroup(A);
   GG:=ConvertGroupToNumberField(G);
   detGG:=sub<GL(1,K)|[GL(1,K)![Determinant(x)]: x in Generators(GG)]>;
   index:=#detGG;
   if index eq 1 then
    A`SLStabilizer:=G;
    return G;
   else
    H:=sub<GG|[x: x in Generators(GG) | Determinant(x) eq 1]>;
    while not #H eq #G/index do
     g:=Random(GG);
     if Determinant(g) eq 1 and not g in H then
      H:=sub<GG|Generators(H) join {g}>;
     end if;
    end while;
    A`SLStabilizer:=ConvertGroupToIntegers(H);
    return A`SLStabilizer;
   end if;
    
  end if;
 else
  G:=HermitianAutomorphismGroup(A);
  GG:=ConvertGroupToNumberField(G);
  order:=#{x: x in GG | CheckMembership(x)};
  H:=sub<GG|[x: x in Generators(GG) | CheckMembership(x)]>;
  while not #H eq order do 
   g:=Random(GG);
   if not g in H and CheckMembership(g) then
    H:=sub<GG|Generators(H) join {g}>;
   end if;
  end while;
  return ConvertGroupToIntegers(H);
 end if;
end intrinsic;

//The following are internal methods, therefore written as functions

matbas:=function(GG);
 //Internal method
 //Convert (2n)*(2n) matrice over Z into n*n matrices over O_K
 //Input&Output: List
 if #GG eq 0 then return []; end if;
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

intrinsic ConvertGroupToNumberField(G::Grp) -> Grp
 {Converts a 2*n by 2*n matrix group over the integers into the corresponding n by n matrix group over K}
 Generators:=SetToIndexedSet(Generators(G));
 ZGENS:=[Generators[i] : i in [1..#Generators]];
 OKGENS:=matbas(ZGENS);
 OG:=sub<GL(n,K)|OKGENS>;
 return OG;
end intrinsic;

intrinsic ConvertGroupToIntegers(G::Grp) -> Grp
 {Converts an n by n matrix group over K into the corresponding 2*n by 2*n matrix group over the integers.}
 Generators:=SetToIndexedSet(Generators(G));
 OKGENS:=[Generators[i] : i in [1..#Generators]];
 ZGENS:=matbas2(OKGENS);
 ZG:=sub<GL(2*n,Integers())|ZGENS>;
 return ZG;
end intrinsic;

intrinsic IsInGL(A::Mtrx) -> BoolElt
 {Checks the input matrix A for membership in GL(L), L as specified in the BasicData of the package.}
 //tests whether A is in GL(L)
 bool:=true;
 for i in [1..n-1] do
  for j in [1..n-1] do
   if not (A[i,j] in Integers(K)) then
    bool:=false;
   end if;
  end for;
  if not (A[i,n] in p2) then
   bool:=false;
  end if;
 end for;
 for j in [1..n-1] do
  if not (A[n,j]) in (ideal<Integers(K)|1>/p2) then
   bool:=false;
  end if;
 end for;
 if not A[n,n] in Integers(K) then
  bool:=false;
 end if;
 if not 1/Determinant(A) in Integers(K) then
  bool:=false;
 end if;
 return bool;
end intrinsic;

