import "../../BasicData.m": n,d,steinitz;
import "../../Initialize.m": K,w,tau,p1,p2,ZB,B,BasHermNorm,mmax,F,sqrtd,Injec;
import "../Functions/Groups.m": matbas,matbas2,matbas3;

limit:=3;

/**********************************************************/

OrbitReps:=function(G,L);
a:=#L;
OO:=[**];
O:=[];
S:=[];

b:=0;
for i in [1..#L] do 
  elt:=L[i];
  neu:=true;
  for oo in OO do 
   if elt in oo then 
    neu:=false; 
    break; 
   end if; 
  end for;
  if neu then 
   Append(~O,i);
   o:=Orbit(G,L[i]);
   s:=Stabiliser(G,L[i]);
   
   // t = OrbitRepresentatives von s auf L 
   OOs:=[**];
   Os:=[];
   bs:=0;
   for ts in [1..#L] do
    neus:=true;
    elts:=L[ts];
    for oos in OOs do 
     if elts in oos then
      neus:=false; 
      break; 
     end if; 
    end for;
    if neus then
     Append(~Os,ts);
     os:=Orbit(s,L[ts]);
     Append(~OOs,os); 
     bs +:= #os;
     if bs eq a then 
      break; 
     end if;
    end if;
   end for; // ts
   Append(~S,Os);
   Append(~OO,o); 
   b +:= #o;
   if b eq a then return O,S; end if;
  end if;
end for;
print [a,b];
return O,S;
end function;


/**********************************************************/

intrinsic CanonicalFormOfMinimalClass(A::HermForm) -> HermForm
 {}
 S:=MinimalVectors(A);
 return NewHermForm(&+[HermitianTranspose(s)*s : s in S]);
end intrinsic;

intrinsic AreEquivalentMinimalClasses(A::HermForm,B::HermForm) -> Any
 {}
 AA:=CanonicalFormOfMinimalClass(A);
 BB:=CanonicalFormOfMinimalClass(B);
 return TestIsometry(NewHermForm(AA`Matrix^-1),NewHermForm(BB`Matrix^-1));
end intrinsic;

intrinsic IsWellRounded(A::HermForm) -> Any
 {}
 S:=MinimalVectors(A);
 V:=KMatrixSpace(K,1,n);
 U:=sub<V|[V!x : x in S]>;
 return Dimension(U) eq n;
end intrinsic;

intrinsic MinimalClasses(V::VData) 
 {Computes a system of representatives of the minimal classes in the Voronoi-tesselation corresponding to V.}

 faceneu:=V`faceneu;


perfectlist:=V`PerfectList;
facelist:=V`FacesList;
facevectList:=V`PerpendicularList;
CriticalValueList:=V`CriticalValueList;


FaceTupleList:=[**];                           //This will contain the tuples of faces of codim>1
Representatives:=[*[p: p in perfectlist] *];            //Representatives of all minimal classes
FaceTupleListOfRepresentatives:=[* [] *];      //We need this to check the dimension of the intersection

print "Starting the computation of minimal classes. Please be patient.";

//n=2

if n eq 2 then

 //Generate the tuples
 for i in [1..#perfectlist] do
  FaceTuples:=[**];
  S:=MinimalVectors(perfectlist[i]);
  for j in [1..#facelist[i]] do
   for k in [j+1..#facelist[i]] do
    Intersection:=facelist[i][j] meet facelist[i][k];
    if #Intersection gt 1 then
     if #FindAllPerp([HermitianTranspose(S[l])*S[l] : l in Intersection]) eq 2 then
      Append(~FaceTuples,[j,k]);
      if k gt #facelist[i] then print "Error."; end if;
     end if;
    end if;
   end for;
  end for;
  Append(~FaceTupleList,FaceTuples);
 end for;
 FaceTupleList:= [* FaceTupleList *];       //This is a bit clumsy; now FTL[1] is the data for codim 2

 //Compute corank 1 classes:

 TempList:=[**];
 for i in [1..#perfectlist] do
  for j in [1..#CriticalValueList[i]] do
   Append(~TempList , perfectlist[i]`Matrix+(CriticalValueList[i][j]/2)*facevectList[i][j] );
  end for;
 end for;

 MinClassReps:=[TempList[1]];
 for x in TempList do
  isi:=false;
  for y in MinClassReps do
   if AreEquivalentMinimalClasses(x,y) then
    isi:=true;
   end if;
  end for;
  if not isi then
   Append(~MinClassReps,x);
  end if;
 end for;
  //Test: append herm forms
 Append(~Representatives,[NewHermForm(x): x in MinClassReps]);
 print "Corank 1 done.";

 //Compute classes of corank >= 2

 codim:=2;

 while n^2-codim ge n do
  TempList:=[**];
  for i in [1..#perfectlist] do
   for j in [1..#FaceTupleList[codim-1][i]] do
    T:=MatrixRing(K,n)!perfectlist[i]`Matrix;
    for k in FaceTupleList[codim-1][i][j] do
     T:=T+(CriticalValueList[i][k]/(2*codim))*facevectList[i][k];
    end for;
    Append(~TempList,T);
   end for;
  end for;

  MinClassReps:=[TempList[1]];
  for x in TempList do
   isi:=false;
   for y in MinClassReps do
    if AreEquivalentMinimalClasses(x,y) then
     isi:=true;
    end if;
   end for;
   if not isi then
    Append(~MinClassReps,x);
   end if;
  end for;
  //Test: Apppend herm forms
  Append(~Representatives,[NewHermForm(x): x in MinClassReps]);
  codim:=codim+1;
 end while;
 print "Data assembled.";
 V`FaceTupleList:=FaceTupleList;
 V`MinimalClasses:=Representatives;
end if;


//n=3

if n eq 3 then
 //Generate the tuples
 print "n=3";

 //Compute corank 1 classes:
 print "Starting the computation of codim 1 classes.";
 TempList:=[**];
 for i in [1..#perfectlist] do
  for j in [1..#CriticalValueList[i]] do
   X:=NewHermForm( perfectlist[i]`Matrix+(CriticalValueList[i][j]/2)*facevectList[i][j] );
   if IsWellRounded(X) then
    Append(~TempList , <X,i,[j]> );
   end if;
  end for;
 end for;

 print "TempList done. Starting equivalence testing for codim 1.";

 MinClassReps:=[TempList[1]];
 for x in TempList do
  isi:=false;
  for y in MinClassReps do
   if AreEquivalentMinimalClasses(x[1],y[1]) then
    isi:=true;
   end if;
  end for;
  if not isi then
   Append(~MinClassReps,x);
  end if;
 end for;

 Append(~Representatives,MinClassReps);
 print "Corank 1 done.";
 
 
 codim:=2;
 while n^2-codim ge limit do
  TempList:=[**];
  print "Now doing it for Codim " cat IntegerToString(codim);
  for rep in Representatives[#Representatives] do
   i:=rep[2];
   subface:=rep[3];
   S:=MinimalVectors(perfectlist[i]);
   for j in [1..#facelist[i]] do
    Intersection:=facelist[i][j] meet &meet[facelist[i][k]: k in subface];
    if #Intersection gt 1 and #FindAllPerp([HermitianTranspose(S[l])*S[l] : l in Intersection]) eq codim then 
     newsubface:=Append(subface,j);
     T:=perfectlist[i]`Matrix;
     for k in newsubface do
      T:=T+(CriticalValueList[i][k]/(2*codim))*facevectList[i][k];
     end for;
     T:=NewHermForm(T);
     if IsWellRounded(T) then
      Append(~TempList,<T,i,newsubface>);
     end if;
    end if;
   end for;
  end for;
  print "TempList done. Starting equivalence testing for codim ", codim,". We have ", #TempList, " candidates." ;

  MinClassReps:=[TempList[1]];
  for x in TempList do
   isi:=false;
   for y in MinClassReps do
    if AreEquivalentMinimalClasses(x[1],y[1]) then
     isi:=true;
    end if;
   end for;
   if not isi then
    Append(~MinClassReps,x);
   end if;
  end for;
  Append(~Representatives,MinClassReps);
  print "Corank ", codim, " done.";
  codim+:=1;   
 end while;   
   
  
  
  
  
  
 print "Representatives assembled, computing translation spaces.";
 //Now we compute bases for the translation spaces of the minimal classes since we already have the data.
 for i in [2..#Representatives] do
  for j in [1..#Representatives[i]] do
   F:=Representatives[i][j];
   TSpaceBas:=[(V`PerfectList[F[2]])`Matrix-(V`PerfectNeighbourList[F[2]][k])`Matrix: k in F[3]];
   TSpaceRat:=[Eltseq(y): y in matbas3(TSpaceBas)];
   assert Rank(Matrix(TSpaceRat)) eq i-1;
   F[1]`TSpaceRat:=KSpaceWithBasis(Matrix(TSpaceRat));
   F[1]`TSpaceBasis:=TSpaceBas;
   Representatives[i][j]:=F;
  end for;
 end for;
 V`FaceTupleList:=Representatives;
 V`MinimalClasses:=[Representatives[1]] cat [[x[1]: x in Representatives[j]]: j in [2..#Representatives]];
end if;
end intrinsic;
