import "../ImQuad/Functions/Groups.m": matbas;


intrinsic ComputeComplexGL(File::MonStgElt,V::VData)
 {This function will assemble all available information about the well-rounded retract and write it to a file (in GAP-readable format)} 
 if not assigned(V`MinimalClasses) then
  MinimalClasses(V);
 end if;
 FormREPS:=V`MinimalClasses;
 
 print "Computing Well-Rounded-Complex for full GL";
 print "Now computing stabilizers";
 Stabs:=AssembleStabilizers(FormREPS,true,IsInGL);
 print "Now computing even stabilizers";
 Evenstabs:=AssembleEvenStabilizers(FormREPS,true,IsInGL,V);
 print "Now computing boundaries";
 BC:=AssembleBoundaryComponents(FormREPS,IsInGL);
 print "Now computing remaining data";
 Dims:=AssembleDimensions(FormREPS);
 Elts:=AssembleElements(BC,Stabs);
 print "Data computed";
 str:="";
 str cat:= DimensionsToGapString(Dims);
 str cat:= ElementsToGapString(Elts);
 str cat:= StabilizersToGapString(Stabs,false);
 str cat:= StabilizersToGapString(Evenstabs,true);
 str cat:= BoundaryComponentsToGapString(BC,Elts);
 str cat:= GeneratorsToGapString(V`ZGens);
 Write(File,str);
end intrinsic;	


intrinsic ComputeComplexLowIndexSubgroup(File::MonStgElt,Reps::Any,CheckMembership::Any,V::VData)
 {This function will assemble all available information about the well-rounded retract for a subgroup of finite index in GL given by its index and a method to check for membership and write it to a file (in GAP-readable format)}
 print "Computing Well-Rounded-Complex for subgroup of low index in GL";
 print "Now computing representatives of minimal classes modulo subgroup action";
 if not assigned(V`MinimalClasses) then
  MinimalClasses(V);
 end if;
 FormREPS:=V`MinimalClasses;
 for i in [2..#FormREPS] do
  //FormREPS[i]:=[NewHermForm(x): x in FormREPS[i]];
  FormREPS[i]:=[x: x in FormREPS[i]];

 end for;
 index:=#Reps;
 Forms:=[[]:i in [1..#FormREPS]];
 for i in [1..#FormREPS] do
  for j in [1..#FormREPS[i]] do
   print i,j;
   S:=StabilizerOfMinimalClass(FormREPS[i][j]);
   LS:=LowIndexStabilizer(FormREPS[i][j],CheckMembership);
   r:=Order(S)/Order(LS);
   s:=index/r;
   z:=1;
   ZZ:=New(HermForm);
   ZZ`Matrix:=Reps[z]*FormREPS[i][j]`Matrix*HermitianTranspose(Reps[z]);
   list:=[ZZ];
   
   while z lt index do 
    z:=z+1;
    if {CheckMembership(Reps[y]*s*Reps[z]^(-1)) : y in [1..z-1],s in S} eq {false} then	
     ZZ:=New(HermForm);
     ZZ`Matrix:=Reps[z]*FormREPS[i][j]`Matrix*HermitianTranspose(Reps[z]);
     Append(~list,ZZ);
    end if;
   end while;
   Forms[i]:=Forms[i] cat list;
  end for;
 end for;
 print "Now computing stabilizers";
 Stabs:=AssembleStabilizers(Forms,false,CheckMembership);
 Evenstabs:=AssembleEvenStabilizers(Forms,false,CheckMembership,V);
 print "Now computing boundaries";
 BC:=AssembleBoundaryComponents(Forms,CheckMembership);
 print "Now computing remaining data";
 Dims:=AssembleDimensions(Forms);
 Elts:=AssembleElements(BC,Stabs);
 print "Data computed"; 
 str:="";
 str cat:= DimensionsToGapString(Dims);
 str cat:= ElementsToGapString(Elts);
 str cat:= StabilizersToGapString(Stabs,false);
 str cat:= StabilizersToGapString(Evenstabs,true);
 str cat:= BoundaryComponentsToGapString(BC,Elts);
 str cat:= GeneratorsToGapString(V`ZGens);
 Write(File,str); 
end intrinsic;



intrinsic ComputeComplexLowIndexSubgroupFixtest(File::MonStgElt,Reps::Any,CheckMembership::Any,V::VData)
 {This function will assemble all available information about the well-rounded retract for a subgroup of finite index in GL given by its index and a method to check for membership and write it to a file (in GAP-readable format)}
 print "Computing Well-Rounded-Complex for subgroup of low index in GL";
 print "Now computing representatives of minimal classes modulo subgroup action";
 if not assigned(V`MinimalClasses) then
  MinimalClasses(V);
 end if;
 FormREPS:=V`MinimalClasses;
 for i in [2..#FormREPS] do
  //FormREPS[i]:=[NewHermForm(x): x in FormREPS[i]];
  FormREPS[i]:=[x: x in FormREPS[i]];

 end for;
 index:=#Reps;
 Forms:=[[]:i in [1..#FormREPS]];
 for i in [1..#FormREPS] do
  for j in [1..#FormREPS[i]] do
   print i,j;
   S:=StabilizerOfMinimalClass(FormREPS[i][j]);
   LS:=LowIndexStabilizer(FormREPS[i][j],CheckMembership);
   r:=Order(S)/Order(LS);
   s:=index/r;
   z:=1;
   ZZ:=New(HermForm);
   ZZ`Matrix:=Reps[z]*FormREPS[i][j]`Matrix*HermitianTranspose(Reps[z]);
   list:=[ZZ];
   
   while z lt index do 
    z:=z+1;
    if {CheckMembership(Reps[y]*s*Reps[z]^(-1)) : y in [1..z-1],s in S} eq {false} then	
     ZZ:=New(HermForm);
     ZZ`Matrix:=Reps[z]*FormREPS[i][j]`Matrix*HermitianTranspose(Reps[z]);
     Append(~list,ZZ);
    end if;
   end while;
   Forms[i]:=Forms[i] cat list;
  end for;
 end for;
 print [#x: x in Forms];
 print "Now computing stabilizers";
 Stabs:=AssembleStabilizers(Forms,false,CheckMembership);
 Evenstabs:=AssembleEvenStabilizers(Forms,false,CheckMembership,V);
 print "Now computing boundaries";
 BC:=AssembleBoundaryComponents(Forms,CheckMembership);
 print "Now computing remaining data";
 Dims:=AssembleDimensions(Forms);
 Elts:=AssembleElements(BC,Stabs);
 print forall{x : x in Elts | CheckMembership(matbas([x])[1])};
 print "Data computed"; 
 str:="";
 str cat:= DimensionsToGapString(Dims);
 str cat:= ElementsToGapString(Elts);
 str cat:= StabilizersToGapString(Stabs,false);
 str cat:= StabilizersToGapString(Evenstabs,true);
 str cat:= BoundaryComponentsToGapString(BC,Elts);
 str cat:= GeneratorsToGapString(V`ZGens);
 Write(File,str); 
end intrinsic;













