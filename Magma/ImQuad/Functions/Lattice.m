import "../../BasicData.m": n,d,steinitz;
import "../../Initialize.m": K,w,tau,p1,p2,ZB,B,BasHermNorm,mmax; 
import "Groups.m": matbas,matbas2,matbas3;

intrinsic TraceForm(A::HermForm) -> Mtrx
 {For a HermForm A returns the trace form of A}
 if assigned A`traceform then return A`traceform; end if;
 res:=MatrixRing(Rationals(),2*n) ! 0;
 for i in [1..2*n] do
  for j in [1..2*n] do
   res[i][j]:=Rationals()!Trace(K!(B[i]*A`Matrix*HermitianTranspose(B[j]))[1][1]);
  end for;
 end for;
 res:=(1/2)*res;
 A`traceform:=res;
 return res;
end intrinsic;

intrinsic EvaluateVector(A::HermForm,x::Mtrx) -> FldQuad
 {Evaluates the vector x at a HermForm A.}
 x:=KMatrixSpace(K,1,n)!x;
 z:=K!0; N:=K!0; I:=ideal<Integers(K)|0>;
 z:=K!((x*(A`Matrix)*HermitianTranspose(x))[1][1]);
 for i in [1..n-1] do
  I:=I+x[1][i]*p1;
 end for;
 I:=I+x[1][n]*(p1/p2);
 N:=Norm(I);
 return z/N;
end intrinsic;

intrinsic ShortenVectors(A::HermForm,m::FldQuadElt,S::SeqEnum) -> SeqEnum
 {for a list of short vectors in a Z-lattice of dimension 2*n by 2*n 
  return the corresponding list of shortest vectors of A with minimum m in dimension n}
 res:=[];
 for i in [1..#S] do
  x:=Vector(S[i][1]);
  xk:=KMatrixSpace(K,1,n)!0;
  for j in [1..2*n] do
   xk:=xk+x[j]*B[j];
  end for;
  if EvaluateVector(A,xk) eq m then
   Append(~res,xk);
  end if;
 end for;
 return res;
end intrinsic;

intrinsic ShortenVector(x::Mtrx) -> Mtrx
 {shortens a 2*n-vector into an n-vector}
 res:=KMatrixSpace(K,1,n)!0;
 for j in [1..2*n] do
  res:=res+x[j]*B[j];
 end for;
 return res;
end intrinsic;

intrinsic IdealNorm(x::Mtrx) -> FldRatElt
 {return the norm of the ideal corresponding to the vector x} 
 N:=0;
 I:=ideal<Integers(K) | 0>;
 for i in [1..n-1] do
  I:=I+x[1][i]*(p1);
 end for;
 I:=I+x[1][n]*(p1/p2);
 N:=Norm(I);
 return N;
end intrinsic;

intrinsic HermitianMinimum(A::HermForm) -> FldRatElt
 {Returns the minimum of the HermForm A on the lattice L.}
 if assigned A`Minimum then return A`Minimum; end if;
 L:=LatticeWithGram(TraceForm(A));
 minL:=Minimum(L);
 S:=ShortVectors(L,minL/mmax,minL*mmax);
 m:=Min([EvaluateVector(A,ShortenVector(s[1])) : s in S]);
 A`Minimum:=m;
 return m;
end intrinsic;

intrinsic PerfectionRankList(M::SeqEnum) -> RngIntElt
 {returns the perfection rank of the set of vectors M}
 VV:=[];
 for m in M do s:=Matrix(m[1]);
  v:=HermitianTranspose(s)*Matrix(s);
  Append(~VV, ElementToSequence(v));
 end for;
 return Rank(Matrix(VV)) ;
end intrinsic;

intrinsic PerfectionRankListOfMatrices(M::SeqEnum) -> RngIntElt
 {returns the perfection rank of the set of matrices M}
 VV:=[];
 for m in M do
  Append(~VV, ElementToSequence(m));
 end for;
 return Rank(Matrix(VV)) ;
end intrinsic;

intrinsic PerfectionRank(M::HermForm) -> RngIntElt
 {Returns the perfection rank of the HermForm M}
 if assigned M`PerfectionRank then return M`PerfectionRank; end if;
 L:=LatticeWithGram(TraceForm(M));
 m:=HermitianMinimum(M);
 S:=ShortVectors(L,m/mmax,m*mmax);
 Sk:=ShortenVectors(M,m,S);
 r:=PerfectionRankList(Sk);
 M`PerfectionRank:=r;
 return r;
end intrinsic;

intrinsic PerfectionCorank(M::HermForm) -> RngIntElt
 {Returns the perfection corank of the Herm Form M.}
 return n^2-PerfectionRank(M);
end intrinsic;

RemoveMultiples:=function(M);
 V:=VectorSpace(K,n);
 out:=[];
 Append(~out,M[1]);
 for m in M do;
 ismultiple:=false;
  for v in out do;
   if Vector(m) in sub<V|[Vector(v)]> then;
    ismultiple:=true;
   end if;
  end for;
  if not ismultiple then;
   Append(~out,m);
  end if;
 end for;
 return out;
end function;

intrinsic MinimalVectors(A::HermForm) -> SeqEnum
 {Returns a mutiple-free list of minimal vectors of the HermForm A.}
 if assigned A`MinimalVectors then return A`MinimalVectors; end if;
 L:=LatticeWithGram(TraceForm(A));
 m:=HermitianMinimum(A);
 S:=ShortVectors(L,m/mmax,m*mmax);
 Sk:=ShortenVectors(A,m,S);
 Sk:=RemoveMultiples(Sk);
 A`MinimalVectors:=Sk;
 return Sk;
end intrinsic;

intrinsic TestIsometry(M::HermForm,N::HermForm : SL:=false , CheckMembership:=0) -> BoolElt,Any
 {Tests A and B for isometry, the result (if there is an isometry) is a matrix g (in rational representation) such that gMg^dagger=N, if SL is set then we test for isometry in SL}
 if CheckMembership cmpeq 0 then
  Me:=TraceForm(M);
  Ne:=TraceForm(N);
  mul1:=Lcm([Denominator(x) : x in Eltseq(Me)]);
  mul2:=Lcm([Denominator(x) : x in Eltseq(Ne)]);
  Me:=ChangeRing(mul1*mul2*Me,Integers());
  Ne:=ChangeRing(mul1*mul2*Ne,Integers());
  LM:=LatticeWithGram(Me);
  LN:=LatticeWithGram(Ne);
  
  Me2:=TraceForm(NewHermForm(tau*M`Matrix));
  Ne2:=TraceForm(NewHermForm(tau*N`Matrix)); 
 
  mul1:=Lcm([Denominator(x) : x in Eltseq(Me2)]);
  mul2:=Lcm([Denominator(x) : x in Eltseq(Ne2)]);
  Me2:=ChangeRing(mul1*Me2,Integers());
  Ne2:=ChangeRing(mul2*Ne2,Integers());
 
  a,b:=IsIsometric(LM,[Me2],LN,[Ne2]);
  
  if not a then
   return false,_;
  elif not SL then 
   return a,b;
  else
   bb:=matbas([b])[1];
   detb:=Determinant(bb);
   if detb eq 1 then
    return true,b;
   end if;
   GG:=ConvertGroupToNumberField(HermitianAutomorphismGroup(M));
   detGG:=sub<GL(1,K)|[GL(1,K)![Determinant(x)]: x in Generators(GG)]>;
   if not GL(1,K)![detb] in detGG then
    return false,_;
   else 
    for x in GG do 
     if Determinant(bb*x) eq 1 then
      return true,matbas2([bb*x])[1];
     end if;
    end for;
   end if;
  end if;
 else
  a,b:=TestIsometry(M,N);
  if not a then 
   return false,_;
  else
   bb:=matbas([b])[1];
   GG:=ConvertGroupToNumberField(HermitianAutomorphismGroup(M));
   for x in GG do
    if CheckMembership(bb*x) then
     return true,matbas2([bb*x])[1];
    end if;
   end for;
   return false,_;
  end if;
 end if;
end intrinsic;


