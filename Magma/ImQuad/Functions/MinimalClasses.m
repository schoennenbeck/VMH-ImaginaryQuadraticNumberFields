import "../../BasicData.m": n,d,steinitz;
import "../../Initialize.m": K,w,tau,p1,p2,ZB,B,BasHermNorm,mmax,F,sqrtd;

intrinsic CanonicalFormOfMinimalClass(F::AlgMatElt) -> AlgMatElt
 {Returns the canonical form T_C for the minimal class C of F.}
 T:=MatrixRing(K,n)!0;
 M:=MinimalVectors(NewHermForm(F));
 for x in M do
  T:=T+HermitianTranspose(x)*x;
 end for;
 return T;
end intrinsic;

intrinsic CanonicalFormOfMinimalClass(M::SeqEnum)-> AlgMatElt
 {Returns the canonical Form of the minimal class given by the list of shortest Vectors}
 T:=MatrixRing(K,n)!0;
 for x in M do
  T:=T+HermitianTranspose(x)*x;
 end for;
 return T;
end intrinsic;

intrinsic StabilizerOfMinimalClass(F::HermForm) -> Grp
 {Returns the automorphism group of the minimal class of F.}
 if assigned F`MinClassStabilizer then
  return F`MinClassStabilizer;
 else
  F`MinClassStabilizer:=ConvertGroupToNumberField(HermitianAutomorphismGroup(NewHermForm(CanonicalFormOfMinimalClass(F)`Matrix^(-1))));
  return StabilizerOfMinimalClass(F);
 end if;
end intrinsic;

intrinsic StabilizerOfMinimalClass(F::AlgMatElt) -> Grp
 {Returns the automorphism group of the minimal class of F.}
 return ConvertGroupToNumberField(HermitianAutomorphismGroup(NewHermForm(CanonicalFormOfMinimalClass(F)^(-1))));
end intrinsic;

intrinsic StabilizerOfMinimalClass(M::SeqEnum) -> Grp
 {Returns the automorphism group of the minimal class given by a list of shortest vectors.}
 return ConvertGroupToNumberField(HermitianAutomorphismGroup(NewHermForm(CanonicalFormOfMinimalClass(M)^(-1))));
end intrinsic;


intrinsic AreEquivalentMinimalClasses(A::AlgMatElt,B::AlgMatElt) -> BoolElt
 {Tests whether the two classes represented by A and B are in the same orbit of the action of GL(L).}
 return TestIsometry(NewHermForm(CanonicalFormOfMinimalClass(A)^(-1)),NewHermForm((CanonicalFormOfMinimalClass(B))^(-1)));
end intrinsic;

intrinsic AreEquivalentMinimalClasses(M::SeqEnum, N::SeqEnum) -> BoolElt
 {Tests whether the two classes given by the minimal vectors M,N are in the same orbit of the action of GL(L).}
 return TestIsometry(NewHermForm((CanonicalFormOfMinimalClass(M))^(-1)),NewHermForm((CanonicalFormOfMinimalClass(N))^(-1)));
end intrinsic;

intrinsic ElementsOfNorm(norm::RngIntElt,P::RngQuadIdl) -> SeqEnum
 {A procedure to produce all elements of given norm in an ideal P over the integers of a number field.}
 OK:=Integers(K);
 Gram:=MatrixRing(Integers(),2)![[2,Trace(tau)],[Trace(tau),2*Norm(tau)]];
 L:=LatticeWithGram(Gram);
 S:=ShortVectors(L,2*norm,2*norm);
 output:=[s[1][1]+s[1][2]*tau: s in S];
 output:=[x: x in output | x in P];
 return output;
end intrinsic;

intrinsic AllMinimalVectors(F::AlgMatElt) -> SeqEnum
 {Computes "all", i.e. including multiples, minimal vectors of the input form F.}
 S:=MinimalVectors(F);
 output:={s: s in S};
 for v in S do
  v10:=true;
  if v[1][1] ne 0 then
   X:=[w/v[1][1]: w in ElementsOfNorm(Norm(v[1][1]),p1)];
   v10:=false;
  end if;
  if v10 then
   X:=[w/v[1][2]: w in ElementsOfNorm(Norm(v[1][2]),p2)];
  else
   X:=[x:x in X| x*v[1][2] in p2];
  end if;
  output:=output join {x*v: x in X};
 end for;
 output:=output join {-x : x in output};
 return [v: v in output];
end intrinsic;

intrinsic AllMinimalVectorsList(LL::SeqEnum) -> SeqEnum
 {Computes "all", i.e. including multiples, minimal vectors from a list of vectors.}
 S:=LL;
 output:={s: s in S};
 for v in S do
  v10:=true;
  if v[1][1] ne 0 then
   X:=[w/v[1][1]: w in ElementsOfNorm(Norm(v[1][1]),p1)];
   v10:=false;
  end if;
  if v10 then
   X:=[w/v[1][2]: w in ElementsOfNorm(Norm(v[1][2]),p2)];
  else
   X:=[x:x in X| x*v[1][2] in p2];
  end if;
  output:=output join {x*v: x in X};
 end for;
 output:=output join {-x : x in output};
 return [v: v in output];
end intrinsic;

intrinsic ReynoldsProjection(G::Grp,A::AlgMatElt) -> AlgMatElt
 {Returns the value of A under the Reynolds operator of G.}
 res:=MatrixRing(K,n)!0;
 for g in G do
  res:=res+g*A*HermitianTranspose(g);
 end for;
 return (1/(#G))*res;
end intrinsic;

intrinsic IsWellRoundedFace(L::SeqEnum) -> BoolElt
 {Takes as input a list of vectors and returns true if the list contains a K-basis of K^n.}
 V:=KMatrixSpace(K,1,n);
 U:=sub<V|L>;
 return Dimension(U) eq n;
end intrinsic;

intrinsic IsWellRoundedForm(A::AlgMatElt) -> BoolElt
 {Checks if the form given by the input matrix A is well-rounded.}
 return IsWellRoundedFace(MinimalVectors(A));
end intrinsic;

