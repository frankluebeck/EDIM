# test.g file for the EDIM package
# `gap4 -b -q -r <test.g' should produce only `true' lines (after the first)
#


ReadPackage("edim", "tst/mat2");
Print("From here you should only see a sequence of `true' lines,\nuntil ",
      "the session quits.\n");

LoadPackage("edim");
#SetInfoLevel(InfoEDIM, 1);


inv := InverseRatMat(mat2);;
inv*mat2=IdentityMat(Length(mat2));
ExponentSquareIntMatFullRank(mat2)=15120;
ElementaryDivisorsPPartRk(mat2, 2)=[ 34, 13, 10, 10, 0 ];
ElementaryDivisorsPPartRkI(mat2, 2, 34)=[ 34, 13, 10, 10, 0 ];
ElementaryDivisorsPPartRkII(mat2, 2, 34)=[ 34, 13, 10, 10, 0 ];
ElementaryDivisorsPPartRkExp(mat2, 2, 34, 5)=[ 34, 13, 10, 10, 0 ];
ElementaryDivisorsPPartRkExpSmall(mat2, 2, 34, 5, 0)=[ 34, 13, 10, 10, 0 ];
ElementaryDivisorsPPartRkExpSmall(mat2+(2^70+1)*2^6,2,34,5,0)=[34,13,10,10,0];
ElementaryDivisorsPPartRkExpSmall(mat2-(2^70+1)*2^6,2,34,5,0)=[34,13,10,10,0];
ElementaryDivisorsSquareIntMatFullRank(mat2)=eldiv2;
ElementaryDivisorsIntMatDeterminant(mat2, 
        1406938906943787632903941535325707304960)=eldiv2;
ElementaryDivisorsPPartHavasSterling(mat2, 2, 5)=[ 34, 13, 10, 10 ];
GcdexIntLLL(21314,345345,564564,768678,42424,64564647,-1313)=
  [ 1, [ -5, -13, -2, 8, -10, 0, 0 ] ];
tr:=HermiteIntMatLLLTrans(mat2);;
tr[2]*mat2=tr[1];
HermiteIntMatLLL(mat2)=tr[1];
tr:=SmithIntMatLLLTrans(mat2);;
tr[2]*mat2*tr[3]=tr[1];
SmithIntMatLLL(mat2)=tr[1];

Add(mat2,Sum(mat2{[1..10]}));
Add(mat2,0*mat2[1]+1);

ElementaryDivisorsPPartRk(mat2, 2)=[ 33, 12, 9, 9, 0 ];
ElementaryDivisorsPPartRkI(mat2, 2, 34)=[ 33, 12, 9, 9, 0 ];
ElementaryDivisorsPPartRkII(mat2, 2, 34)=[ 33, 12, 9, 9, 0 ];
ElementaryDivisorsPPartRkExp(mat2, 2, 34, 5)=[ 33, 12, 9, 9, 0 ];
ElementaryDivisorsPPartRkExpSmall(mat2, 2, 34, 5, 0)=[ 33, 12, 9, 9, 0 ];
ElementaryDivisorsPPartHavasSterling(mat2, 2, 5)=[ 33, 12, 9, 9 ];
tr:=HermiteIntMatLLLTrans(mat2);;
tr[2]*mat2=tr[1];
HermiteIntMatLLL(mat2)=tr[1];
tr:=SmithIntMatLLLTrans(mat2);;
tr[2]*mat2*tr[3]=tr[1];
SmithIntMatLLL(mat2)=tr[1];

mat2 := MutableTransposedMat(mat2);;

ElementaryDivisorsPPartRk(mat2, 2)=[ 33, 12, 9, 9, 0 ];
ElementaryDivisorsPPartRkI(mat2, 2, 34)=[ 33, 12, 9, 9, 0 ];
ElementaryDivisorsPPartRkII(mat2, 2, 34)=[ 33, 12, 9, 9, 0 ];
ElementaryDivisorsPPartRkExp(mat2, 2, 34, 5)=[ 33, 12, 9, 9, 0 ];
ElementaryDivisorsPPartRkExpSmall(mat2, 2, 34, 5, 0)=[ 33, 12, 9, 9, 0 ];
ElementaryDivisorsPPartHavasSterling(mat2, 2, 5)=[ 33, 12, 9, 9 ];
tr:=HermiteIntMatLLLTrans(mat2);;
tr[2]*mat2=tr[1];
HermiteIntMatLLL(mat2)=tr[1];
tr:=SmithIntMatLLLTrans(mat2);;
tr[2]*mat2*tr[3]=tr[1];
SmithIntMatLLL(mat2)=tr[1];
quit;
