###########################################################################
##  
#F  RankMod( <A>, <p> ) . . . . . . . rank of integer matrix <A> modulo <p>
##  
##  Here <p> must not be a prime.
##  
DeclareGlobalFunction("RankMod2");
InstallGlobalFunction(RankMod2, function(arg)
  local   A,  p,  A1, Tp, A1Tp, n,  m, p2,  inv,  res,  first,  v,  vv,  vp,
          i,  c,  pos,  x, y, pr, j, l, gcdrep, gcd, row, irows, submatrix;
  
  A := arg[1];
  p := arg[2];
  if Length(arg) > 2 then
    submatrix := true;
  fi;
  
  # use fast 8bit vectors if p small
  if p<256  and IsPrimeInt(p) and not submatrix then
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
  Tp := [1..n];
  A1Tp := [];
  inv := [];
  res := [];
  irows := [];

  for row in [1..Length(A)] do
    v := A[row];
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
      Add(irows, row);
    fi;
  od;
  if submatrix then
    return  [Length(A1Tp), irows, Tp];
  else
    return Length(A1Tp);
  fi;
end);
