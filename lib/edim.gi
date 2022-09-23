#############################################################################
##
#A  edim.gi                  EDIM package                        Frank Lübeck
##
##  
#Y  Copyright (C) 1999  Lehrstuhl D f\"ur Mathematik, RWTH Aachen
##  
##  This file contains the  main functions of  the EDIM package.  The main
##  reference for the algorithms implemented here is
##  
##  Frank Lübeck, On the  Computation  of Elementary Divisors  of  Integer
##  Matrices,  to appear in Journal of Symbolic Computation.
##  
##  A preprint version is available from
##  https://www.math.rwth-aachen.de/~Frank.Luebeck/preprints/preprints_en.html
##  
  
###########################################################################
##  
#F  PAdicLinComb(mat, invp, p, r, v) . . . . . . . .  
#F  PAdicLinComb([res, mat, invp, p, rr, vv], newr ). . . . . . . . p-adic
#F  approximation of solution of mat*x=v
##  
##  mat must be an integer matrix, p a prime,  invp an inverse modulo p of
##  mat, v a vector on which mat operates and  r some natural number. This
##  functions   returns a  list of six   objects  [res, mat, invp, p,  rr,
##  vv]. Here res is the (unique) integer vector with entries in the range
##  ]-p^r/2..p^r/2] such that  vv = mat*res has  only entries divisible by
##  p^r. mat, invp and p are as in the input.  rr=r if vv<>0. In case vv=0
##  rr is the  smallest number  such that the  entries  of res are in  the
##  range ]-p^rr/2..p^rr/2].
##  
##  In the  second  form the first  argument  is the output  of a previous
##  call.  The result  is the same  as  that of PAdicLinComb(mat, invp, p,
##  newr, v) if newr>=rr.
##  
InstallGlobalFunction(PAdicLinComb, function(arg)
  local   mat, invp, p, r0, r, v, res,  p2,  i,  x,  j,  y, pi;
  if Length(arg)>2 then
    mat := arg[1]; invp := arg[2]; p := arg[3]; 
    r := arg[4]; v := -arg[5]; 
    res := 0*[1..Length(mat)];
    r0 := 0;
  else
    x := arg[1];
    mat := x[2]; invp := x[3]; p := x[4]; 
    r := arg[2]; v := x[6]; r0 := x[5];
    res := x[1];
  fi;
  
  p2 := QuoInt(p,2);
  pi := p^r0;
  for i in [r0..r-1] do
    if v=0*v then
      return [res, mat, invp, p, i+1, v];
    fi;
    x := -v*invp;
    for j in [1..Length(x)] do
      y := x[j] mod p;
      if y>p2 then
        x[j] := y-p;
      else
        x[j] := y;
      fi;
    od;
    res := res + pi*x;
    pi := pi*p;
    v := (x*mat+v)/p;
  od;
  return [res, mat, invp, p, r, v];
end);


###########################################################################
##  
#F  InverseRatMat( mat[, p ] ) . . . .  . . .  returns inverse of a matrix 
#F  over the rationals
##  
##  This function computes  the inverse of an  invertible matrix over  the
##  rationals using p-adic approximations.
##    
##  This  seems  to be better  than  the standard Gauss algorithm (mat^-1)
##  already          for     small         matrices.           (Try,  e.g,
##  RandomMat(20,20,[-10000..10000]) or RandomMat(100,100).)  The optional
##  argument p must be a prime such that mat mod  p is invertible (default
##  is p=251).
##  
InstallGlobalFunction(InverseRatMat,  function(arg)
  local   m,  p,  n,  nr,  den,  mip,  index,  null,  inv,  i,  v,  
          nf,  vv,  rr,  prr,  prr2,  c,  d,  j;
  m := arg[1];
  if Length(arg)>1 then
    p := arg[2];
  else
    p := 251;
  fi;
  n := Length(m);
  
  den := Lcm(Set(List(Concatenation(m),DenominatorRat)));
  if den <> 1 then
    m := m*den;
  fi;
  mip := fail;
  while mip=fail do
    Info(InfoEDIM, 1, "Inverting matrix modulo ",p,".\n");
    mip := InverseIntMatMod(m, p);
    if mip=fail then p := NextPrimeInt(p); fi;
  od;
  index := 1;
  null := 0*[1..n];
  inv := [1..n];
  for i in [1..n] do
    Info(InfoEDIM, 1, i," \c");
    v := ShallowCopy(null);
    v[i] := index;
    nf := true;
    vv := [ShallowCopy(null), m, mip, p, 0, -v];
    rr := LogInt(index,p)+2;
    prr := p^rr;
    prr2 := p^QuoInt(rr,2);
    while nf do
      vv := PAdicLinComb(vv, rr);
      if vv[6]=null then
        Info(InfoEDIM, 1, "+\c");
        nf := false;
        inv[i] := vv[1]/index;  
      else
##          c := List(vv[1], x-> RatNumber(prr, QuoInt(prr2,2), prr2, x));
##  this line is substituted by more complicated ones which save a lot
##  of computations if the approximation is not yet sufficient.         
        c := [RatNumberFromModular(prr, QuoInt(prr2,2), prr2, vv[1][1])];
        if IsRat(c[1]) then
          d := DenominatorRat(c[1]);
        fi;
        j := 1;
        while IsRat(c[j]) and j<n do
          j := j+1;
          c[j] := RatNumberFromModular(prr, QuoInt(prr2,2), prr2, (vv[1][j]*d)
                          mod prr);
          if c[j]<>fail then
            c[j] := c[j]/d;
            d := LcmInt(d, DenominatorRat(c[j]));
          fi;
        od;
        if not IsVector(c) then
          Info(InfoEDIM, 1, "-\c");
          rr := rr + 10;
          prr := p^rr;
          prr2 := p^QuoInt(rr,2);
        else
          d := Lcm(List(c, DenominatorRat));
          if  (c*d)*m = d*v then
            nf := false;
            Info(InfoEDIM, 1, " (",d,") ");
            inv[i] := c/index;
            index := index * d;
          else
            Info(InfoEDIM, 1, "#\c");
            rr := rr + 10;
            prr := p^rr;
            prr2 := p^QuoInt(rr,2);
          fi;
        fi;
      fi;
    od;
    Info(InfoEDIM, 1, "[depth: ",vv[5],"] ");
  od;
  Info(InfoEDIM, 1, "\n");
  if den <> 1 then
    return inv*den;
  else
    return inv;
  fi;
end);

###########################################################################
##  
#F  RationalSolutionIntMat( mat, v[, p[, imat ] ] ) . . . . returns a
#F  solution x with rational entries of x * mat = v
##  
##  mat must be an invertible (over Q) matrix with integer entries and
##  v a vector of integers of the same length as mat. The optional arguments
##  are a prime p such that mat mod p is invertible and a matrix imat over
##  GF(p) which is the inverse of mat mod p.
##  
##  The solution is constructed by p-adic approximation, as in 
##  InverseRatMat.
##  
##  args: m, v[, p, mip] 
InstallGlobalFunction(RationalSolutionIntMat,  function(arg)
  local   m,  v,  p,  mip,  n,  null,  nf,  vv,  rr,  prr,  prr2,  c,  
          d,  j;
  m := arg[1]; v := arg[2];
  if Length(arg)>2 then  
    p := arg[3]; 
  else 
    p := 251; 
  fi;
  if Length(arg)>3 then
    mip := arg[4];    
  else
    mip := InverseIntMatMod(m,  p); 
  fi; 
  
  n := Length(m);
  null := 0*[1..n];
  nf := true;
  vv := [0*[1..n], m, mip, p, 0, -v];
  rr := 2*LogInt(Maximum(List(v, AbsInt))+1, p) + 2;
  prr := p^rr;
  prr2 := p^QuoInt(rr,2);
  while nf do
    vv := PAdicLinComb(vv, rr);
    if vv[6]=null then
      Info(InfoEDIM, 1, "+\c");
      nf := false;
      return [vv[1], 1];  
    else
      ##      c := List(vv[1], x-> RatNumber(prr, QuoInt(prr2,2), prr2, x));
      ##  this line is substituted by more complicated ones which save a lot
      ##  of computations if the approximation is not yet sufficient.         
      c := [RatNumberFromModular(prr, QuoInt(prr2,2), prr2, vv[1][1])];
      if IsRat(c[1]) then
        d := DenominatorRat(c[1]);
      fi;
      j := 1;
      while IsRat(c[j]) and j<n do
        j := j+1;
        c[j] := RatNumberFromModular(prr, QuoInt(prr2,2), prr2, (vv[1][j]*d)
                        mod prr);
        if c[j]<>fail then
          c[j] := c[j]/d;
          d := LcmInt(d, DenominatorRat(c[j]));
        fi;
      od;
      if not IsVector(c) then
        Info(InfoEDIM, 1, "-\c");
        rr := rr + 10;
        prr := p^rr;
        prr2 := p^QuoInt(rr,2);
      else
        d := Lcm(List(c, DenominatorRat));
        if  (c*d)*m = d*v then
          nf := false;
          Info(InfoEDIM, 1, " (",d,") ");
          return [c, d];
        else
          Info(InfoEDIM, 1, "#\c");
          rr := rr + 10;
          prr := p^rr;
          prr2 := p^QuoInt(rr,2);
        fi;
      fi;
    fi;
  od;
end);

###########################################################################
##  
#F  ExponentSquareIntMatFullRank( mat[, p[, nr ] ] ) . . . returns biggest 
#F  elementary divisor
##  
##  For  a  square integer matrix <mat>   of  full rank   the least common
##  multiple of all entries of the inverse  matrix <mat>^-1 is exactly the
##  biggest elementary divisor of <mat>.
##    
##  This is just  a   slight modification of 'InverseRatMat'.    The third
##  argument nr tells the function to return the lcm  of the first nr rows
##  of the rational  inverse matrix  only. Very  often the function   will
##  already  return the biggest elementary  divisor  with  nr=2 or 3  (and
##  spend some time in checking, that this is correct).
##  
##  This function  runs into an error if  the matrix is not invertible mod
##  p.
##  
InstallGlobalFunction(ExponentSquareIntMatFullRank, function(arg)
  local   m, p, nr, mip,  bnd,  r,  index,  null,  i,  n, v,  nf,  vv,  
          rr,  prr,  prr2,  c,  d, j;
  m := arg[1];
  if Length(arg)>1 then
    p := arg[2];
  else
    p := 251;
  fi;
  n := Length(m);
  if Length(arg)>2 then
    nr := arg[3];
  else
    nr := n;
  fi;
  
  mip := fail;
  while mip=fail do
    Info(InfoEDIM, 1, "Inverting matrix modulo ",p,".\n");
    mip := InverseIntMatMod(m, p);
    if mip=fail then
      p := NextPrimeInt(p);
    fi;
  od;
  
  index := 1;
  null := 0*[1..n];
  for i in [1..nr] do
    Info(InfoEDIM, 1, i," \c");
    v := ShallowCopy(null);
    v[i] := index;
    nf := true;
    vv := [ShallowCopy(null), m, mip, p, 0, -v];
    rr := LogInt(index,p)+2;
    prr := p^rr;
    prr2 := p^QuoInt(rr,2);
    while nf do
      vv := PAdicLinComb(vv, rr);
      if vv[6]=null then
        Info(InfoEDIM, 1, "+\c");
        nf := false;
      else
        c := [RatNumberFromModular(prr, QuoInt(prr2,2), prr2, vv[1][1])];
        if IsRat(c[1]) then
          d := DenominatorRat(c[1]);
        fi;
        j := 1;
        while IsRat(c[j]) and j<n do
          j := j+1;
          c[j] := RatNumberFromModular(prr, QuoInt(prr2,2), prr2, (vv[1][j]*d)
                          mod prr);
          if c[j]<>fail then
            c[j] := c[j]/d;
            d := LcmInt(d, DenominatorRat(c[j]));
          fi;
        od;
        if not IsVector(c) then
          Info(InfoEDIM, 1, "-\c");
          rr := rr + 10;
          prr := p^rr;
          prr2 := p^QuoInt(rr,2);
        else
          d := Lcm(List(c, DenominatorRat));
          if  (c*d)*m = d*v then
            nf := false;
            Info(InfoEDIM, 1, " (",d,") ");
            index := index * d;
          else
            Info(InfoEDIM, 1, "#\c");
            rr := rr + 10;
            prr := p^rr;
            prr2 := p^QuoInt(rr,2);
          fi;
        fi;
      fi;
    od;
    Info(InfoEDIM, 1, "[depth: ",rr,"] ");
  od;
  Info(InfoEDIM, 1, "\n");
  return index;
end);

###########################################################################
##  
##  The following  programs are the main part  of this  mini-package. They
##  implement an algorithm which finds for a  given prime p the p-parts of
##  the elementary  divisors of an integer  matrix.  As  input the rank of
##  the matrix is needed.  If available a bound for the highest power of p
##  dividing the biggest elementary divisor can be given.
##  
##  The algorithm is explained in:
##  
##  Frank L\"ubeck, On the Computation  of Elementary Divisors of  Integer
##  Matrices, submitted. (Available as preprint-file upon request from the
##  author.)
##  

##########################################################################
##  
#F  ElementaryDivisorsPPartRk( <A>, <p>[, <rk> ])  . . . . . . . . .
#F  ElementaryDivisorsPPartRkI( <A>, <p>, <rk> )  . . . . . . . . .
#F  ElementaryDivisorsPPartRkII( <A>, <p>, <rk> )  . . . . . . . . .
#F  ElementaryDivisorsPPartRkExp( <A>, <p>, <rk>, <exp> ) . . . . . . 
#F  ElementaryDivisorsPPartRkExpSmall( <A>, <p>, <rk>, <exp>, <il> ) . . . 
#F  returns list [m_1, m_2,  .., m_r]  where m_i is the number of nonzero 
#F  elementary divisors of A divisible by p^i
##  
##  
##  <A> must be a  matrix with integer entries,  <p> a prime, and <rk> the
##  rank of <A> (as rational matrix).
##  
##  In the form with  <exp>  must be an  upper bound for the highest power
##  of <p>  appearing in an  elementary  divisor of <A>. This  information
##  allows  reduction  of   matrix entries  modulo   <p>^<exp> during  the
##  computation.
##  
##  Here <exp> is  allowed  to be too   small. In this case  the  function
##  returns 'fail'.
##  
##  The  command  with  the  'Small'  ending   can  be  used  as  long  as
##  <p>^(<exp>+1) is  an immediate  integer and  (<p>^(<exp>-1))(p-1) fits
##  into an unsigned C long integer.  This is a kernel function, which can
##  be fast even for large matrices.
##  
InstallGlobalFunction(ElementaryDivisorsPPartRk, function(arg)
  local   A,  p,  m,  n,  rk,  r,  res, tmp, z;
  A := arg[1];
  p := arg[2];
  m := Length(A);
  n := Length(A[1]);
  if Length(arg)=2 then
    # rank not given, we first check mod 251 if full rank
    Info(InfoEDIM, 1, "Compute rank mod 251 . . .\n");
    tmp := A*Z(251)^0;
    for z in tmp do CONV_VEC8BIT(z, 251); od;
    rk := RankMat(tmp);
    if not rk=Minimum(n, m) then
      Info(InfoEDIM, 1, "Compute rank in characteristic 0 . . .\n");      
      rk := RankMat(A);
    fi;
  else
    rk := arg[3];
  fi;
  r := Maximum(3, LogInt(2^31-1, p)-2);
  res := fail;
  while res=fail do
    Info(InfoEDIM, 1, "Calling ElementaryDivisorsPPartRkExp with exp ",r," . . .\n");
    res := ElementaryDivisorsPPartRkExp(A, p, rk, r);
    r := r+3;
  od;
  return res;
end);

InstallGlobalFunction(ElementaryDivisorsPPartRkI, function(A, p, rk)
  local   A1,  A2,  A2l, Tp, pr,  m,  n,  i,  j,  inv,  res,  i0,  ii,  
          x,  c,  vv,  pos,  i1,  r, max;
  
  m := Length(A);
  n := Length(A[1]);
  Tp := [1..n];
  A2 := List(A, ShallowCopy);
  
  A1 := [];
  inv := [];
  res := [];

  while Length(A1)<rk  do
    i0 := Length(A2);
    A2l := 0;
    for ii in [1..i0] do
      vv := A2[ii];
      Unbind(A2[ii]);
      for i in [1..Length(A1)] do
	c := (vv[Tp[i]]*inv[i]) mod p;
        if c <> 0 then
          vv := vv - c*A1[i];
	fi;
      od;
      pos := Length(A1) + 1;
      while pos<=n and vv[Tp[pos]] mod p = 0 do
        pos := pos+1;
      od;
      if pos>n then
	Info(InfoEDIM, 1, "-\c");
        A2l := A2l+1;
        A2[A2l] := vv/p;  
      else
	Info(InfoEDIM, 1, "+\c");
        i1 := Length(A1)+1;
        A1[i1] := vv;
 
	if pos <> i1 then
          x := Tp[pos];
          Tp[pos] := Tp[i1];
          Tp[i1] := x;
        fi;
	inv[i1] := (1/(A1[i1][Tp[i1]] mod p)) mod p;
      fi;
    od;
    Add(res, rk-Length(A1));
    Info(InfoEDIM, 1, "\n# max(all entries) = ", Maximum(List(Concatenation(
         Concatenation(A1),Concatenation(A2)), AbsInt)), "\n#Rank found: ",
         Length(A1), "\n");
  od;
  return res;
end);

## comment on this
InstallGlobalFunction(ElementaryDivisorsPPartRkII,  function(A, p, rk)
  local   A1,  A1m, A2,  A2n, Tp, m,  n,  p2,  inv,  res,  ll,  i0,  nrnew,  
          lcnew,  ii,  vv, lc,  c, i,  nr,  x,  len,  v,  pos,  i1,  max;
  
  m := Length(A);
  n := Length(A[1]);
  Tp := [1..n];
  A2 := List(A, ShallowCopy);
  p2 := QuoInt(p,2);
  
  A1 := [];
  A1m := [];
  inv := [];
  res := [];

  while Length(A1)<rk  do
    ll := Length(A1);
    i0 := Length(A2);
    A2n := [];
    nrnew := [];
    lcnew := [];
    for ii in [1..i0] do
      vv := List(A2[ii], x-> x mod p);
      lc := 0*[1..Length(A1)];
      for i in [1..Length(A1)] do
	c := (vv[Tp[i]]*inv[i]) mod p;
        if c <> 0 then
          vv := vv - c*A1[i];
          if i>ll then
            nr := nrnew[i-ll];
            x := lcnew[nr];
            len := Length(x);
            lc{[1..len]} := lc{[1..len]}-c*x;
          else
            lc[i] := -c;
          fi;
        fi;
      od;
      vv := List(vv, x-> x mod p);
      for i in [1..Length(lc)] do
        x := lc[i] mod p;
        if x>p2 then
          lc[i] := x-p;
        else
          lc[i] := x;
        fi;
      od;
      if ll>0 then
        v := A2[ii] + lc{[1..ll]}*A1{[1..ll]};
      else
        v := A2[ii];
      fi;
      if Length(lc)>ll then
        v := v + lc{[ll+1..Length(lc)]}*A2{nrnew};
      fi;
      pos := Length(A1) + 1;
      while pos<=n and vv[Tp[pos]] = 0 do
        pos := pos+1;
      od;
      if pos>n then
	Info(InfoEDIM, 1, "-\c");
        Add(A2n, v/p);
      else
	Info(InfoEDIM, 1, "+\c");
        i1 := Length(A1)+1;
        A1m[i1] := vv;
        A1[i1] := v;
        Add(lc, 1);
        lcnew[ii] := lc;
        Add(nrnew,ii);
	if pos <> i1 then
          x := Tp[pos];
          Tp[pos] := Tp[i1];
          Tp[i1] := x;
        fi;
	inv[i1] := (1/(A1[i1][Tp[i1]] mod p)) mod p;
      fi;
    od;
    A2 := A2n;
    Add(res, rk-Length(A1));
    Info(InfoEDIM, 1, "\n# max(all entries) = ", Maximum(List(Concatenation(
         Concatenation(A1),Concatenation(A2)), AbsInt)), "\n#Rank found: ",
         Length(A1), "\n");
  od;
  return res;
end);

InstallGlobalFunction(ElementaryDivisorsPPartRkExp,  function(A, p, rk, r)
  local   A1,  A2,  A2l, Tp, pr,  b, m,  n,  i,  j,  inv,  res,  i0,  ii,  
          x,  c,  vv,  pos,  i1;
  
  pr := p^(r+1);
  
  # if everything small then delegate to kernel function
  b := 8 * GAPInfo.BytesPerVariable;
  if IsBoundGlobal("ElementaryDivisorsPPartRkExpSmall") and (p-1)*(pr-1) < 2^b 
     and IsPlistRep(A) and ForAll(A, IsPlistRep) and pr < 2^(b-4) then
    return ValueGlobal("ElementaryDivisorsPPartRkExpSmall") 
           (A, p, rk, r, InfoLevel(InfoEDIM));
  fi;
  m := Length(A);
  n := Length(A[1]);
  Tp := [1..n];
  A2 := List(A, ShallowCopy);
  
  A1 := [];
  inv := [];
  res := [];

  while Length(A1)<rk and r>=0 do
    i0 := Length(A2);
    A2l := 0;
    for ii in [1..i0] do
      vv := A2[ii];
      Unbind(A2[ii]);
      for i in [1..Length(A1)] do
	c := (vv[Tp[i]]*inv[i]) mod p;
        if c <> 0 then
          AddRowVector(vv, A1[i], -c);
	fi;
      od;
      for j in [1..n] do
        vv[j] := vv[j] mod pr;
      od;
      pos := Length(A1) + 1;
      while pos<=n and vv[Tp[pos]] mod p = 0 do
        pos := pos+1;
      od;
      if pos>n then
	Info(InfoEDIM, 1, "-\c");
        A2l := A2l+1;
        A2[A2l] := vv/p;  
      else
	Info(InfoEDIM, 1, "+\c");
        i1 := Length(A1)+1;
        A1[i1] := vv;
 
	if pos <> i1 then
          x := Tp[pos];
          Tp[pos] := Tp[i1];
          Tp[i1] := x;
        fi;
	inv[i1] := (1/A1[i1][Tp[i1]]) mod p;
      fi;
    od;
    Add(res, rk-Length(A1));
    Info(InfoEDIM, 1, "\n#Rank found: ", Length(A1), "\n");
    pr := pr/p;
    r := r-1;
  od;
  if Length(A1)=rk then
    return res;
  else
    return fail;
  fi;
end);

###########################################################################
##  
##  Two functions to put things together.
##    
#F  ElementaryDivisorsSquareIntMatFullRank( mat )  . . elementary divisors
##    
##  Here   first  the   biggest  elementary  divisor     is computed   via
##  'ExponentSquareIntMatFullRank'. If it  runs into an error  because mat
##  is singular mod a  choosen prime (it starts  with  251) just  type two
##  times 'return;' - then it tries the next prime.
##  
##  The rest is done using 'ElementaryDivisorsPPartRkExp' and 'RankMod'.
##  
##  The  functions   fails if  the  biggest  elementary divisor  cannot be
##  completely factored and the non-factored part is  not a divisor of the
##  biggest elementary divisor only.
##  
##  
#F  ElementaryDivisorsIntMatDeterminant( mat, det[, rk] ) . . . elementary 
#F  divisors where the biggest determinant divisor is given
##  
##  det can be given as integer in the form of Collected(FactorsInt(det)).
##  If the matrix does not have full rank then its rank must be given.
##  
InstallGlobalFunction(ElementaryDivisorsSquareIntMatFullRank, function(arg)
  local   mat, index,  fac,  n,  res,  i,  a,  rk,  mul;
  
  mat := arg[1];
  Info(InfoEDIM, 1, "Computing largest elementary divisor . . .\n");
  index := CallFuncList(ExponentSquareIntMatFullRank,arg);
  
  Info(InfoEDIM, 1, "(Partly) factoring the largest elementary divisor . . .\n");
  fac := Collected(CheapFactorsInt(index));
  n := Length(mat);  
  
  # to make it GAP4 usable
  res := 0*[1..n]+1;
  for a in fac do
    if IsPrimeInt(a[1]) then
      if a[2]=1 then
        Info(InfoEDIM, 1, "For ",a[1],"-part of elementary divisors only rank must ",
              "be computed . . .\c");
        rk := RankMod(mat, a[1]);
        Info(InfoEDIM, 1, "(corank ", n-rk, ")\n");
        for i in [rk+1..n] do 
          res[i] := res[i]*a[1];
        od;
      else
        Info(InfoEDIM, 1, "For ",a[1],"-part of elementary divisors the modular ",
              "algorithm is used . . .\n");
        mul := ElementaryDivisorsPPartRkExp(mat, a[1], n, a[2]);
        for rk in n-mul do
          for i in [rk+1..n] do
            res[i] := res[i]*a[1];
          od;
        od;
      fi;
    else 
      Info(InfoEDIM, 1, "Checking if non-factored part is only in last elementary",
            " divisor . . .\n");
      rk := RankMod(mat, a[1]);
      if rk = n-1 or (IsList(rk) and ForAll(rk, a-> a[2] = n-1)) then
        Info(InfoEDIM, 1, "Success.\n");
        res[n] := res[n]*a[1]^a[2];
      else
        Info(InfoEDIM, 1, "Not successful, going back to library function.\n");
        return ElementaryDivisorsMat(mat);
      fi;
    fi;
  od;
  return res;
end);

InstallGlobalFunction(ElementaryDivisorsIntMatDeterminant,  function(arg)
  local   mat,  det,  n,  fac,  res,  i,  a,  rk,  mul;

  mat := arg[1];
  det := arg[2];
  if Length(arg)>2 then
    n := arg[3];
  else
    n := Minimum(Length(mat), Length(mat[1]));
  fi;
  
  if IsInt(det) then
    Info(InfoEDIM, 1, "(Partly) factoring the determinant . . .\n");
    fac := Collected(CheapFactorsInt(det));
  else
    fac := det;
  fi;
  
  res := 0*[1..n]+1;
  for a in fac do
    if IsPrimeInt(a[1]) then
      if a[2] = 1 then
        res[n] := res[n]*a[1];
      elif a[2]<=3 then
        # only three possibilities which can be distinguished by rank mod p 
        Info(InfoEDIM, 1, "For ",a[1],"-part of elementary divisors only rank must ",
              "be computed . . .\n");
        rk := RankMod(mat, a[1]);
        for i in [rk+1..n] do 
          res[i] := res[i]*a[1];
        od;
        res[n] := res[n]*a[1]^(a[2]-n+rk);
      else
        Info(InfoEDIM, 1, "For ",a[1],"-part of elementary divisors the modular ",
              "algorithm is used . . .\n");
        mul := ElementaryDivisorsPPartRkExp(mat, a[1], n, a[2]);
        for rk in n-mul do
          for i in [rk+1..n] do
            res[i] := res[i]*a[1];
          od;
        od;
      fi;
    else 
      Info(InfoEDIM, 1, "Checking if non-factored part is only in last elementary",
            " divisor . . .\n");
      rk := RankMod(mat, a[1]);
      if rk = n-1 or (IsList(rk) and ForAll(rk, a-> a[2] = n-1)) then
        Info(InfoEDIM, 1, "Success.\n");
        res[n] := res[n]*a[1]^a[2];
      else
        Info(InfoEDIM, 1, "Not successful, going back to library function.\n");
        return ElementaryDivisorsMat(mat);
      fi;
    fi;
  od;
  return res;
end);

###########################################################################
##  
##  The following function  implements the modular algorithm  described in
##  
##  Havas, Sterling: Integer matrices and abelian groups, Springer LNCS 72 
##
##  We added a slight improvement: we divide the considered submatrices by
##  the p-part of the  gcd of all  entries  (and lower the d  below). This
##  reduces the size of the entries and often shortens the pivot search.
##  
#F  ElementaryDivisorsPPartHavasSterling( A, p, d ) . . .  . . . . returns
#F  the list of p-parts of the nonzero elementary divisors of A
##  
##  Here a lower bound  d for the highest  power of p dividing the biggest
##  elementary divisor     of A must  be  given.   Smaller  d  improve the
##  performance of the algorithm considerably.
##  
InstallGlobalFunction(ElementaryDivisorsPPartHavasSterling, function(mat, p, d)
  local   pd,  m,  n,  r,  inv,  div,  min,  i,  j,  k,  x,  pk,  u, res;
  # we need a power of p higher than that in largest elementary divisor
  d := d+1;
  pd := p^(d);
  m := Length(mat);
  n := Length(mat[1]);
  r := 1;
  mat := List(mat, ShallowCopy);
  # for collecting the p-parts of the elementary divisors
  inv := [];
  div := 1;
  
  while true do
    # we consider submatrix [r..m]x[r..n]
    Info(InfoEDIM, 1, r, " \c");
    
    # find entry with smallest multiplicity of p, stop when an entry
    # prime to p was found
    min := [0,0,d+1];
    i := r;
    j := r;
    while min[3]>0 and i<=m and j <= n do
      if mat[i][j] <> 0 then
        k := 0;
        x := mat[i][j]/p;
        while k < d and IsInt(x) do 
          k := k+1;
          x := x/p;
        od;
        if k = d then
          mat[i][j] := 0;
        fi;
        if k < min[3] then
          min := [i,j,k];
        fi;
      fi;
      if j = n then
        j := 1;
        i := i+1;
      else
        j := j+1;
      fi;
    od;
    if min[1] = 0 then
      # ready if all entries are divisible by p^d
      res := [];
      inv := Filtered(inv, a-> a<>1);
      while Length(inv)<>0 do
        Add(res, Length(inv));
        inv := Filtered(inv/p, a-> a<>1);
      od;
      return res; 
    fi;
    # swap rows and columns if necessary
    if min[1] <> r then
      k := min[1];
      x := mat[r];
      mat[r] := mat[k];
      mat[k] := x;
    fi;
    if min[2] <> r then
      k := min[2];
      for i in [r..m] do
        x := mat[i][r];
        mat[i][r] := mat[i][k];
        mat[i][k] := x;
      od;
    fi;
    
    # divide everything by common power of p
    if min[3] > 0 then
      pk := p^min[3];
      div := div * pk;
      d := d - min[3];
      pd := p^d;
      for i in [r..m] do
        for j in [r..n] do
          mat[i][j] := mat[i][j]/pk;
        od;
      od;
    fi;
    Add(inv, div);
    Info(InfoEDIM, 1, "(",div,")\c");
    
    # make mat[r][r] = 1 mod p^d
    u := Gcdex(mat[r][r], pd).coeff1;
    if AbsInt(2*u) > pd then
      u := u - SignInt(u)*pd;
    fi;
    for j in [r..n] do
      mat[r][j] := (u*mat[r][j]) mod pd;
      if 2*mat[r][j] > pd then
        mat[r][j] := mat[r][j] - pd;
      fi;
    od;
    
    # reduce column
    Info(InfoEDIM, 1, ".\c");
    for i in [r+1..m] do
      u := mat[i][r];
      mat[i][r] := 0;
      for j in [r+1..n] do
        mat[i][j] := (mat[i][j] - u*mat[r][j]) mod pd;
        if 2*mat[i][j] > pd then
          mat[i][j] := mat[i][j] - pd;
        fi;
      od;
    od;
    
    # forget r-th row and column
    r := r + 1;
  od;
end);
