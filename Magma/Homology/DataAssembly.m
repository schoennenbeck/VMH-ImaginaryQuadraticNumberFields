import "../BasicData.m": n;
import "../Initialize.m": K,tau;

intrinsic AssembleDimensions(Forms::Any)-> SeqEnum
 {Returns the list of dimensions of the well-rounded complex}
 return [#Forms[i]:i in [1..#Forms]];
end intrinsic;

intrinsic AssembleStabilizers(Forms::Any,gl::BoolElt,CheckMembership::Any)-> SeqEnum
 {Returns the list of stabilizers of cells either in GL or in a LowIndexSubgroup}
 if gl then
  Stabs:=[[Generators(ConvertGroupToIntegers(StabilizerOfMinimalClass(Forms[i][j]))): j in [1..#Forms[i]]]: i in [1..#Forms]];
 else
  Stabs:=[[Generators(ConvertGroupToIntegers(LowIndexStabilizer(Forms[i][j],CheckMembership))): j in [1..#Forms[i]]]: i in [1..#Forms]];
 end if;
 return Stabs;
end intrinsic;

intrinsic AssembleEvenStabilizers(Forms::Any,gl::BoolElt,CheckMembership::Any,V::VData) -> SeqEnum
 {Same as above for the orientation preserving stabilizers}
 if gl then 
  EvenStabs:=[[Generators(ConvertGroupToIntegers(EvenStabilizerOfMinimalClass(Forms[i][j],V))): j in [1..#Forms[i]]]: i in [1..#Forms]];
 else
 EvenStabs:=[[Generators(ConvertGroupToIntegers(LowIndexEvenStabilizer(Forms[i][j],CheckMembership,V))): j in [1..#Forms[i]]]: i in [1..#Forms]];
 end if;
 return EvenStabs;
end intrinsic;

intrinsic AssembleBoundaryComponents(Forms::Any,CheckMembership::Any) -> List
 {Returns the list of Boundaries as needed later}
 BoundaryComponents:=[**];
 for j in [2..#Forms] do
  tmpj:=[**];
  for k in [1..#Forms[j]] do
   tmpk:=[**];
   for i in [1..#Forms[j-1]] do
    tmpi:=[**];
    for g in BoundaryEmbeddings(Forms[j][k],Forms[j-1][i],CheckMembership) do
     Append(~tmpi, [*i,g*]);
    end for;
    if #tmpi ne 0 then 
     tmpk:=tmpk cat tmpi;
    end if;
   end for;
   if #tmpk ne 0 then
    Append(~tmpj,tmpk);
   end if;
  end for;
  if #tmpj ne 0 then
   Append(~BoundaryComponents,tmpj);
  end if;
 end for;
 return BoundaryComponents;
end intrinsic;

intrinsic AssembleElements(BoundaryComponents::List,Stabilizers::SeqEnum)->List
 {Assembles all needed elements}
 B:=BoundaryComponents;
 S:=Stabilizers;
 
 elts:=[**];
 for i in [1..#S] do
  for j in [1..#S[i]] do
   for x in MatrixGroup<2*n,Integers()|S[i][j]> do	
    if Position(elts,x) eq 0  then	
     Append(~elts,x);
    end if;
   end for;
  end for;
 end for;
 for i in [1..#B] do
  for j in [1..#B[i]] do
   for k in [1..#B[i][j]] do
    for x in Generators(ConvertGroupToIntegers(MatrixGroup<n,K|B[i][j][k][2]>)) do
     if Position(elts,x) eq 0 then
      Append(~elts,x);
     end if;
    end for;
   end for;
  end for;
 end for;
 l:=#elts;
 /*for i in [1..l] do
  for j in [1..l] do
   if Position(elts,elts[i]*elts[j]) eq 0 then
    Append(~elts,elts[i]*elts[j]);
   end if;
  end for;
 end for;*/
 return elts;
end intrinsic;
