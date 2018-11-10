import "BasicData.m": n,d,steinitz;

K<w>:=QuadraticField(d);

if d mod 4 eq 1 then
 tau:=(1+w)/2;
else
 tau:=w;
end if;
C,f:=ClassGroup(K);
Idealvertreter:=[];
for c in C do Append(~Idealvertreter,f(c)); end for;

mmax:=Maximum([Norm(p) : p in Idealvertreter]);

if d ne -1 then
 F<sqrtd>:=QuadraticField(-d);
 Injec:=hom<F -> RealField() | Sqrt(-d)>;
else
 F:=Rationals();
 sqrtd:=1;
 Injec:=hom<F -> RealField() |>;
end if;

//Choose a suitable representative for a (the considered lattice is O_K^{n-1} + a) from Idealvertreter
p1:=Idealvertreter[1];
p2:=Idealvertreter[steinitz];

//Z-base for a:
if p2 eq p1 then ZB:=[1,tau]; else ZB:=Basis(p2); end if;

//Z-base for  L

B:=[];
for k in [1..(n-1)] do
 v:=KMatrixSpace(K,1,n)!0;
 for i in [1..2] do
  if i eq 1 then
   v[1][k]:=K!1;
  else
   v[1][k]:=tau;
  end if;
  Append(~B,v);
 end for;
end for; 
 
for i in [1..2] do
 v:=KMatrixSpace(K,1,n)!0;
 if i eq 1 then
  v[1][n]:=ZB[1];
 else
  v[1][n]:=ZB[2];
 end if;
 Append(~B,v);
end for;

//Determine normalized R-base of space of Hermitian forms

BasHermNorm:=[];
for i in [1..n] do
 res:=MatrixRing(K,n)!0;
 res[i][i]:=1;
 Append(~BasHermNorm,res);
 for j in [i+1..n] do
  for k in [1..2] do
   res:=MatrixRing(K,n)!0;
   if k eq 1 then
    res[i][j]:=1/2;
    res[j][i]:=1/2;
    Append(~BasHermNorm,res);
   else
    res[i][j]:=w/2;
    res[j][i]:=-w/2;
    Append(~BasHermNorm,res);
   end if;
  end for;
 end for;
end for;





























