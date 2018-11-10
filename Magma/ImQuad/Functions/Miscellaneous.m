import "../../BasicData.m": n,d,steinitz;
import "../../Initialize.m": K,w,tau,p1,p2,ZB,B,BasHermNorm,mmax,F,sqrtd;

intrinsic SymmetricCoordinates(A::Mtrx) -> Any
 {For the Hermitian matrix A return its coordinates w.r.t. the basis of Hermitian matrices specified in the initialization file.}
 V:=KMatrixSpaceWithBasis([KMatrixSpace(K,n,n)!x : x in BasHermNorm]);
 return ChangeRing(Vector(Coordinates(V,V!A)),Rationals());
end intrinsic;

intrinsic CoordinatesToMatrix(L::Any) -> Mtrx
 {Returns a Hermitian matrix which has the entered coordinates w.r.t. the basis of Hermitian matrices specified in the initialization file.}
 return &+[L[i]*BasHermNorm[i] : i in [1..#BasHermNorm]];
end intrinsic;

intrinsic FindPerp(L::SeqEnum) -> Mtrx
 {Returns a Hermitian matrix which is perpendicular to the Hermitian matrices entered in the input list L.}
 M:=KMatrixSpace(Rationals(),#BasHermNorm,#L)!0;
 for i in [1..#BasHermNorm] do
  for j in [1..#L] do
   M[i][j]:=Rationals()!Trace(BasHermNorm[i]*L[j]);
  end for;
 end for;
 k:=KernelMatrix(M);
 z:=k[1];

 return &+[z[i]*BasHermNorm[i] : i in [1..#BasHermNorm]];
end intrinsic;

intrinsic FindAllPerp(L::SeqEnum) -> Mtrx
 {Returns a Hermitian matrix which is perpendicular to the Hermitian matrices entered in the input list L.}
 M:=KMatrixSpace(Rationals(),#BasHermNorm,#L)!0;
 for i in [1..#BasHermNorm] do
  for j in [1..#L] do
   M[i][j]:=Rationals()!Trace(BasHermNorm[i]*L[j]);
  end for;
 end for;
 k:=KernelMatrix(M);
 
 return [&+[k[j][i]*BasHermNorm[i] : i in [1..#BasHermNorm]]: j in [1..NumberOfRows(k)]];
end intrinsic;



intrinsic RealToString(r::FldReElt) -> MonStgElt
 {}
 if Sign(r) eq -1 then
  str := "-";
 else
  str := "";
 end if;
 r:=Abs(r);
 p := Integers()! Floor(r) ;
 str := str cat IntegerToString(p) cat ".";
 for i := 1 to 15 do
  r:=10*(r-p);
  p := Integers()! Floor(r) ;
  str := str cat IntegerToString(p);
 end for;
 return str;
end intrinsic;

