gap> START_TEST("EDIM");
gap> if not IsBound(ElementaryDivisorsPPartRkExpSmall) then
>    Print("    Call this test file only with properly installed EDIM   !!!\n");
>    Print("    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
> fi;
gap> ReadPackage("edim", "tst/mat2");;
Reading 34x34 integer matrix 'mat2' with elementary divisors 'eldiv2'.
gap> inv := InverseRatMat(mat2);;
gap> inv*mat2=IdentityMat(Length(mat2));
true
gap> ExponentSquareIntMatFullRank(mat2)=15120;
true
gap> ElementaryDivisorsPPartRk(mat2, 2)=[ 34, 13, 10, 10, 0 ];
true
gap> ElementaryDivisorsPPartRkI(mat2, 2, 34)=[ 34, 13, 10, 10, 0 ];
true
gap> ElementaryDivisorsPPartRkII(mat2, 2, 34)=[ 34, 13, 10, 10, 0 ];
true
gap> ElementaryDivisorsPPartRkExp(mat2, 2, 34, 5)=[ 34, 13, 10, 10, 0 ];
true
gap> ElementaryDivisorsPPartRkExpSmall(mat2, 2, 34, 5, 0)=[ 34, 13, 10, 10, 0 ];
true
gap> ElementaryDivisorsSquareIntMatFullRank(mat2)=eldiv2;
true
gap> ElementaryDivisorsIntMatDeterminant(mat2, 
>         1406938906943787632903941535325707304960)=eldiv2;
true
gap> ElementaryDivisorsPPartHavasSterling(mat2, 2, 5)=[ 34, 13, 10, 10 ];
true
gap> GcdexIntLLL(21314,345345,564564,768678,42424,64564647,-1313)=
>   [ 1, [ -5, -13, -2, 8, -10, 0, 0 ] ];
true
gap> tr:=HermiteIntMatLLLTrans(mat2);;
gap> tr[2]*mat2=tr[1];
true
gap> HermiteIntMatLLL(mat2)=tr[1];
true
gap> tr:=SmithIntMatLLLTrans(mat2);;
gap> tr[2]*mat2*tr[3]=tr[1];
true
gap> SmithIntMatLLL(mat2)=tr[1];
true
gap> ElementaryDivisorsPPartRkExpSmall(mat2+(2^70+1)*2^6, 2, 34, 5, 0)=
> [34,13,10,10,0];
true
gap> ElementaryDivisorsPPartRkExpSmall(mat2-(2^70+1)*2^6, 2, 34, 5, 0)=
> [34,13,10,10,0];
true
gap> b:=GAPInfo.BytesPerVariable*8;;
gap> p:=NextPrimeInt(2^(b/2-2)-1000);;
gap> ElementaryDivisorsPPartRkExpSmall(mat2, p, 34, 1, 0);
Error, exponent too large, see ?ElementaryDivisorsPPartRkExpSmall
gap> ElementaryDivisorsPPartRkExpSmall(mat2, 2, 34, b-5, 0);
Error, exponent too large, see ?ElementaryDivisorsPPartRkExpSmall
gap> ElementaryDivisorsPPartRkExpSmall(mat2, 2, 34, b-6, 0);
[ 34, 13, 10, 10, 0 ]
gap> Add(mat2,Sum(mat2{[1..10]}));
gap> Add(mat2,0*mat2[1]+1);
gap> 
gap> ElementaryDivisorsPPartRk(mat2, 2)=[ 33, 12, 9, 9, 0 ];
true
gap> ElementaryDivisorsPPartRkI(mat2, 2, 34)=[ 33, 12, 9, 9, 0 ];
true
gap> ElementaryDivisorsPPartRkII(mat2, 2, 34)=[ 33, 12, 9, 9, 0 ];
true
gap> ElementaryDivisorsPPartRkExp(mat2, 2, 34, 5)=[ 33, 12, 9, 9, 0 ];
true
gap> ElementaryDivisorsPPartRkExpSmall(mat2, 2, 34, 5, 0)=[ 33, 12, 9, 9, 0 ];
true
gap> ElementaryDivisorsPPartHavasSterling(mat2, 2, 5)=[ 33, 12, 9, 9 ];
true
gap> tr:=HermiteIntMatLLLTrans(mat2);;
gap> tr[2]*mat2=tr[1];
true
gap> HermiteIntMatLLL(mat2)=tr[1];
true
gap> tr:=SmithIntMatLLLTrans(mat2);;
gap> tr[2]*mat2*tr[3]=tr[1];
true
gap> SmithIntMatLLL(mat2)=tr[1];
true
gap> mat2 := MutableTransposedMat(mat2);;
gap> 
gap> ElementaryDivisorsPPartRk(mat2, 2)=[ 33, 12, 9, 9, 0 ];
true
gap> ElementaryDivisorsPPartRkI(mat2, 2, 34)=[ 33, 12, 9, 9, 0 ];
true
gap> ElementaryDivisorsPPartRkII(mat2, 2, 34)=[ 33, 12, 9, 9, 0 ];
true
gap> ElementaryDivisorsPPartRkExp(mat2, 2, 34, 5)=[ 33, 12, 9, 9, 0 ];
true
gap> ElementaryDivisorsPPartRkExpSmall(mat2, 2, 34, 5, 0)=[ 33, 12, 9, 9, 0 ];
true
gap> ElementaryDivisorsPPartHavasSterling(mat2, 2, 5)=[ 33, 12, 9, 9 ];
true
gap> tr:=HermiteIntMatLLLTrans(mat2);;
gap> tr[2]*mat2=tr[1];
true
gap> HermiteIntMatLLL(mat2)=tr[1];
true
gap> tr:=SmithIntMatLLLTrans(mat2);;
gap> tr[2]*mat2*tr[3]=tr[1];
true
gap> SmithIntMatLLL(mat2)=tr[1];
true
gap> STOP_TEST("EDIM", 0);
