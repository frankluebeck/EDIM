#############################################################################
##
#A  util.gi              EDIM-mini-package                     Frank Lübeck
##
##  
#Y  Copyright (C) 1999  Lehrstuhl D f\"ur Mathematik, RWTH Aachen
##  
##  In this file we install the utility functions of the EDIM package.
##  

###########################################################################
##  
##  Only slight modification of 'FactorsInt' in library to avoid error
##  message if number not completely factored
##  
#F  CheapFactorsInt( n ) . . . . . 
#F  CheapFactorsInt( n, nr ) . . . . .  . . returns list  of factors of n,
#F  including 'small' prime factors - here nr is  the number of iterations
#F  for 'FactorsRho' (default: 2000)
##  
InstallGlobalFunction(CheapFactorsInt, function ( arg )
    local  n,  nr,  sign,  factors,  p,  tmp;
    
    n := arg[1];
    if Length(arg) > 1 then
      nr := arg[2];
    else
      nr := 2000;
    fi;
    
    # make $n$ positive and handle trivial cases
    sign := 1;
    if n < 0  then sign := -sign;  n := -n;  fi;
    if n < 4  then return [ sign * n ];  fi;
    factors := [];

    # do trial divisions by the primes less than 1000
    # faster than anything fancier because $n$ mod <small int> is very fast
    for p  in Primes  do
        while n mod p = 0  do Add( factors, p );  n := QuoInt(n, p);  od;
        if n < (p+1)^2 and 1 < n  then Add(factors,n);  n := 1;  fi;
        if n = 1  then factors[1] := sign*factors[1];  return factors;  fi;
    od;

    # do trial divisions by known factors
    for p  in Primes2  do
        while n mod p = 0  do Add( factors, p );  n := QuoInt(n, p);  od;
        if n = 1  then factors[1] := sign*factors[1];  return factors;  fi;
    od;

    # handle perfect powers
    p := SmallestRootInt( n );
    if p < n  then
        while 1 < n  do
            Append( factors, CheapFactorsInt(p, nr) );
            n := QuoInt(n, p);
        od;
        Sort( factors );
        factors[1] := sign * factors[1];
        return factors;
    fi;

    # let 'FactorsRho' do the work
    tmp := FactorsRho( n, 1, 16, nr );
    if 0 < Length(tmp[2])  then
        Print( "# sorry,  cannot factor ", tmp[2], "\n" );
    fi;
    factors := Concatenation( factors, tmp[1], tmp[2] );
    Sort( factors );
    factors[1] := sign * factors[1];
    return factors;
end);

###########################################################################
##  
#F  HadamardBoundIntMat( mat ) . . . . . . . . . . . . .  as the name says
##    
##  The Hadamard  bound is the product  of Euclidean  norms of the nonzero
##  rows (or columns) of mat. It is an upper  bound for the absolute value
##  of the determinant of mat.
##  
InstallGlobalFunction(HadamardBoundIntMat, function(m)
  local   s,  a,  x,  b;
  s := 1;
  for a in m do 
    x := a * a;
    if x<>0 then
      s := s * x;
    fi;
  od;
  return RootInt(s);
end);

###########################################################################
##  
#F  RatNumberFromModular(n, k, l, x) . . . . . . . . returns r/s = x mod n
##  
##  n,  k, l must be positive  integers with 2*k*l<=n and  x  should be an
##  integer with -n/2  < x <=  n/2. If it exists  this  function returns a
##  rational number r/s with 0 < s < l, gcd(s, n) = 1, -k  < r < k and r/s
##  congruent to x mod  n (i.e., n| r -s  x).  Such an  r/s is unique. The
##  function returns `fail' if such a number does not exist.
##  
InstallGlobalFunction(RatNumberFromModular, function(n,k,l,x)
  local   s,  a,  b,  as,  bs,  q,  y;
  s := SignInt(x);
  x := s*x;
  a := n;                 b := 0;
  as := x;                bs := 1;
  while as >= k and bs < l do
    q := QuoInt(a,as);    s := -s;
    x := a-q*as;          y := b+q*bs;
    a := as;              b := bs;
    as := x;              bs := y;
  od;
  if as < k and bs <= l then
    return s*as/bs;
  else
    return fail;
  fi;
end);

###########################################################################
##  
#F  InverseIntMatMod( mat, p ) . . . . . . . inverse matrix modulo prime p
##  
##  This functions returns  for an integer matrix  'mat' an integer matrix
##  'inv'  with entries  in  the  range ]-p/2..p/2] such   that inv  * mat
##  reduced modulo p is the identity matrix.
##  
##  It runs into an error if the inverse does not exist.
##  
InstallGlobalFunction(InverseIntMatMod, function(mat, p)
  local   p2,  m,  n,  i,  j,  x;
  p2 := QuoInt(p,2);
  m :=  (mat*Z(p)^0);
  if p=2 then
    for x in m do CONV_GF2VEC(x); od;
  elif p<256 and IsPrimeInt(p) then
    for x in m do CONV_VEC8BIT(x, p); od;
  fi;
  m := m^-1;
  if m = fail then
    return fail;
  fi;
  # XXX temp.hack ! XXX
  m := List(m, ShallowCopy);
  n := Length(mat);
  for i in [1..n] do
    for j in [1..n] do
      x := IntFFE(m[i][j]);
      if x > p2 then
        m[i][j] := x-p;
      else
        m[i][j] := x;
      fi;
    od;
  od;
  return m;
end);

###########################################################################
##  
#F  RankMod( <A>, <p> ) . . . . . . . rank of integer matrix <A> modulo <p>
##  
##  Here <p> must not be a prime.
##  
InstallGlobalFunction(RankMod, function(A, p)
  local   A1, A2, A2n, Tp, A1Tp, n,  m, p2,  inv,  res,  first,  v,  vv,  vp,
          i,  c,  pos,  x, y, pr, j, l, gcdrep, gcd;
  
  # use fast 8bit vectors if p small
  if p<256  and IsPrimeInt(p) then
    Info(InfoEDIM, 1, 1, "Using fast 8bit vectors for `RankMod'.\n");
    m :=  (A*Z(p)^0);
    if p = 2 then
      for x in m do CONV_GF2VEC(x); od;
    else
      for x in m do CONV_VEC8BIT(x, p); od;
    fi;
    return RankMat(m);
  fi;
  
  Info(InfoEDIM, 1, 1, "Using standard algorithm for `RankMod'.\n");  
  n := Length(A[1]);
  p2 := QuoInt(p, 2);
  A2 := A;
  Tp := [1..n];
  A1Tp := [];
  inv := [];
  res := [];

  A2n := [];
  for v in A2 do
    vp := [1..n];
    for j in [1..n] do
      vp[j] := v[Tp[j]] mod p;
    od;
    l := Length(A1Tp); 
    i := 0;
    for i in [1..l] do
      c := (vp[i] * inv[i]) mod p;
      if c>p2 then
        c := c - p;
      fi;

      if c<>0 then
        vp := vp - c * A1Tp[i];
      fi;
    od;
    l := l + 1;
    pos := l;
    while pos<=n and vp[pos] mod p = 0 do
      pos := pos + 1;
    od;
    if pos>n then
      Info(InfoEDIM, 1, "-\c");     
    else
      Info(InfoEDIM, 1, "+\c");     
      for j in [1..n] do
        vp[j] := vp[j] mod p;
      od;
      Add(A1Tp, vp);
      if pos <> l then
        for y in A1Tp do
          x := y[l];
          y[l] := y[pos];
          y[pos] := x;
        od;
        x := Tp[l];
        Tp[l] := Tp[pos];
        Tp[pos] := x;
      fi;
      
      gcdrep := GcdRepresentationOp(Integers, A1Tp[l][l], p);
      gcd := gcdrep * [A1Tp[l][l], p];
      if gcd <> 1 then
        l := Set(Concatenation(CheapFactorsInt(gcd), 
                     CheapFactorsInt(QuoInt(p, gcd))));
        res := List(l, x->RankMod(A, x));
        for i in [1..Length(res)] do
          if IsInt(res[i]) then
            res[i] := [[l[i], res[i]]];
          fi;
        od;
        return Concatenation(res);
      fi;
      Add(inv, gcdrep[1] mod p);
    fi;
  od;
  return Length(A1Tp);  
end);
