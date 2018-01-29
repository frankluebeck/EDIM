#############################################################################
##
#A  hmmlll.gi             EDIM-mini-package                     Frank Lübeck
##
##  
#Y  Copyright (C) 1999  Lehrstuhl D f\"ur Mathematik, RWTH Aachen
##  
##  The following  functions implement    an extended  Gcd-algorithm   for
##  integers and  a  Hermite and Smith  normal  form algorithm for integer
##  matrices using LLL techiques. They are described in the paper
##  
##  Havas,  Majewski,  Matthews: "Extended  gdd  and  Hermite  normal form
##  algorithms  via  lattice basis  reduction", Experimental Mathematics 7
##  (1998) 125-135
##  
##  (The programs don't have  many comments because  these can be found in
##  the clearly written  paper. We apply  the algorithm for Hermite normal
##  form several times  to get the Smith normal  form, that is not in  the
##  paper.)
##  
##  They are particularly useful  if one  wants to  have the  normal forms
##  together with transforming matrices.  These transforming matrices have
##  spectacularly nice (i.e., small) entries.
##  
##  In detail:
##  
#F  GcdexIntLLL(n_1,  n_2, ...) . . .  . . .  . returns for integers n_i a
#F  list [g,  [c_1, c_2, ...]], where g=c_1*n_1+c_2*n_2+...  is the gcd of
#F  the n_i
##  
#F  HermiteIntMatLLL( mat) . . . . . . . . . . returns Hermite normal form
#F  of an integer matrix mat
##  
#F  HermiteIntMatLLLTrans( mat  )  . . . . . . . . . returns  [H, L] where
#F  H=L*mat is the Hermite normal form of an integer matrix mat
##  
#F  SmithIntMatLLL( mat) . . . . . . . . . . . . returns Smith normal form
#F  of an integer matrix mat
##  
#F  SmithIntMatLLLTrans( mat  )  . . . . . . . . . returns [S, L, R] where
#F  S=L*mat*R is the Smith normal form of an integer matrix mat
##  
InstallGlobalFunction(GcdexIntLLL, function(arg)
  local   m,  red1,  swap,  B,  D,  la,  r,  a,  m1,  n1,  k,  i;
  m := Length(arg);
  if m=1 then 
    a := SignInt(arg[1]); 
    if a=0 then a:=1; fi;
    return [a*arg[1], a]; 
  fi;
  
  red1 := function(k, i)
    local   q;
    if a[i] <> 0 then
      q := BestQuoInt(a[k], a[i]);
    elif 2*AbsInt(la[k][i]) > D[i+1] then
      q := BestQuoInt(la[k][i], D[i+1]);
    else
      return;
    fi;
    a[k] := a[k] - q*a[i];
    B[k] := B[k] - q*B[i];
    la[k][i] := la[k][i] - q*D[i+1];
    la[k]{[1..i-1]} := la[k]{[1..i-1]} - q*la[i];
  end;
  
  swap := function(k)
    local   x,  j,  i;
    x := a[k];
    a[k] := a[k-1];
    a[k-1] := x;
    x := B[k];
    B[k] := B[k-1];
    B[k-1] := x;
    for j in [1..k-2] do
      x := la[k][j];
      la[k][j] := la[k-1][j];
      la[k-1][j] := x;
    od;
    for i in [k+1..m] do
      x := la[i][k-1]*D[k+1] - la[i][k]*la[k][k-1];
      la[i][k-1] := (la[i][k-1]*la[k][k-1] + la[i][k]*D[k-1])/D[k];
      la[i][k] := x/D[k];
    od;
    D[k] := (D[k-1]*D[k+1] + la[k][k-1]^2)/D[k];
  end;
  
  
  # initialization
  B := IdentityMat(m);
  la := [1..m];
  for r in [2..m] do
    la[r] := 0*[1..r-1];
  od;
  
  # we shift all indices of D's by one (since we don't have D[0] in GAP)
  D := 1+0*[0..m];
  a := ShallowCopy(arg);
  m1 := 3; n1 := 4;
  k := 2;
  
  # now the LLL routine
  while k <= m do
    red1(k, k-1);
    if a[k-1] <> 0 or (a[k-1] = 0 and a[k] = 0 and
               n1*(D[k-1]*D[k+1]+la[k][k-1]^2) < m1*D[k]^2) then
      swap(k);
      if k > 2 then k := k-1; fi;
    else
      for i in [k-2, k-3..1] do
        red1(k, i);
      od;
      k := k+1;
    fi;
  od;
  if a[m] < 0 then
    a[m] := -a[m];
    B[m] := -B[m];
  fi;
  
  return [a[m], B[m]];
end);


InstallGlobalFunction(HermiteIntMatLLLTrans, function(mat)
  local   m,  n,  red2,  swap,  B,  D,  col1,  col2,  la,  r,  a,  m1,
          n1,  j,  k,  i; 
  m := Length(mat);
  n := Length(mat[1]);
  
  red2 := function(k, i)
    local jj,  q;
    col1 := 1;
    while col1 <= n and a[i][col1] = 0 do
      col1 := col1 + 1;
    od;
    if col1 <= n and a[i][col1] < 0 then
      B[i] := -B[i];
      la[i] := -la[i];
      for jj in [i+1..m] do
        la[jj][i] := -la[jj][i];
      od;
      a[i] := -a[i];
    fi;
    col2 := 1;
    while col2 <= n and a[k][col2] = 0 do
      col2 := col2 + 1;
    od;
      
    if col1 <= n then
      q := BestQuoInt(a[k][col1], a[i][col1]);
    elif 2*AbsInt(la[k][i]) > D[i+1] then
      q := BestQuoInt(la[k][i], D[i+1]);
    else
      return;
    fi;
    a[k] := a[k] - q*a[i];
    B[k] := B[k] - q*B[i];
    la[k][i] := la[k][i] - q*D[i+1];
    la[k]{[1..i-1]} := la[k]{[1..i-1]} - q*la[i];
  end;
  
  swap := function(k)
    local   x,  j,  i;
    x := a[k];
    a[k] := a[k-1];
    a[k-1] := x;
    x := B[k];
    B[k] := B[k-1];
    B[k-1] := x;
    for j in [1..k-2] do
      x := la[k][j];
      la[k][j] := la[k-1][j];
      la[k-1][j] := x;
    od;
    for i in [k+1..m] do
      x := la[i][k-1]*D[k+1] - la[i][k]*la[k][k-1];
      la[i][k-1] := (la[i][k-1]*la[k][k-1] + la[i][k]*D[k-1])/D[k];
      la[i][k] := x/D[k];
    od;
    D[k] := (D[k-1]*D[k+1] + la[k][k-1]^2)/D[k];
  end;
  
  
  # initialization
  B := IdentityMat(m);
  la := [1..m];
  for r in [2..m] do
    la[r] := 0*[1..r-1];
  od;
  
  # we shift all indices of D's
  D := 1+0*[0..m];
  a := List(mat,ShallowCopy);
  m1 := 3; n1 := 4;
  
  # sign adjustment
  j := 1;
  while j<=n and a[m][j]=0  do
    j := j+1;
  od;
  if j <= n and a[m][j] < 0 then
    a[m] := -a[m];
    B[m][m] := -1;
  fi;
  
  k := 2;
  
  # now the LLL routine
  Info(InfoEDIM, 1, "Hermite normal form of matrix with ",m," rows . . .\n");
  while k <= m do
    Info(InfoEDIM, 1, k," \c");
    red2(k, k-1);
    if (col1 <= col2 and col1 <= n)
       or (col1 = col2 and col1 = n+1 and
           n1*(D[k-1]*D[k+1]+la[k][k-1]^2) < m1*D[k]^2) then
      swap(k);
      if k > 2 then k := k-1; fi;
    else
      for i in [k-2, k-3..1] do
        red2(k, i);
      od;
      k := k+1;
    fi;
  od;
  Info(InfoEDIM, 1, "\n");
  # adjust ordering to usual convention
  return [Reversed(a), Reversed(B)];
end);


InstallGlobalFunction(HermiteIntMatLLL, function(mat)
  local   m,  n,  red2,  swap,  D,  col1,  col2,  la,  r,  a,  m1,
          n1,  j,  k,  i; 
  m := Length(mat);
  n := Length(mat[1]);
  
  red2 := function(k, i)
    local jj,  q;
    col1 := 1;
    while col1 <= n and a[i][col1] = 0 do
      col1 := col1 + 1;
    od;
    if col1 <= n and a[i][col1] < 0 then
      la[i] := -la[i];
      for jj in [i+1..m] do
        la[jj][i] := -la[jj][i];
      od;
      a[i] := -a[i];
    fi;
    col2 := 1;
    while col2 <= n and a[k][col2] = 0 do
      col2 := col2 + 1;
    od;
      
    if col1 <= n then
      q := BestQuoInt(a[k][col1], a[i][col1]);
    elif 2*AbsInt(la[k][i]) > D[i+1] then
      q := BestQuoInt(la[k][i], D[i+1]);
    else
      return;
    fi;
    a[k] := a[k] - q*a[i];
    la[k][i] := la[k][i] - q*D[i+1];
    la[k]{[1..i-1]} := la[k]{[1..i-1]} - q*la[i];
  end;
  
  swap := function(k)
    local   x,  j,  i;
    x := a[k];
    a[k] := a[k-1];
    a[k-1] := x;
    for j in [1..k-2] do
      x := la[k][j];
      la[k][j] := la[k-1][j];
      la[k-1][j] := x;
    od;
    for i in [k+1..m] do
      x := la[i][k-1]*D[k+1] - la[i][k]*la[k][k-1];
      la[i][k-1] := (la[i][k-1]*la[k][k-1] + la[i][k]*D[k-1])/D[k];
      la[i][k] := x/D[k];
    od;
    D[k] := (D[k-1]*D[k+1] + la[k][k-1]^2)/D[k];
  end;
  
  
  # initialization
  la := [1..m];
  for r in [2..m] do
    la[r] := 0*[1..r-1];
  od;
  
  # we shift all indices of D's
  D := 1+0*[0..m];
  a := List(mat,ShallowCopy);
  m1 := 101; n1 := 400;
  
  # sign adjustment
  j := 1;
  while j<=n and a[m][j]=0  do
    j := j+1;
  od;
  if j <= n and a[m][j] < 0 then
    a[m] := -a[m];
  fi;
  
  k := 2;
  
  # now the LLL routine
  Info(InfoEDIM, 1, "Hermite normal form of matrix with ",m," rows . . .\n");
  while k <= m do
    Info(InfoEDIM, 1, k," \c");
    red2(k, k-1);
    if (col1 <= col2 and col1 <= n)
       or (col1 = col2 and col1 = n+1 and
           n1*(D[k-1]*D[k+1]+la[k][k-1]^2) < m1*D[k]^2) then
      swap(k);
      if k > 2 then k := k-1; fi;
    else
      for i in [k-2, k-3..1] do
        red2(k, i);
      od;
      k := k+1;
    fi;
  od;
  
  Info(InfoEDIM, 1, "\n");
  # adjust ordering to usual convention
  return Reversed(a);
end);

InstallGlobalFunction(SmithIntMatLLL, function(mat)
  local   m1,  d,  i,  j,  g;
  m1 := HermiteIntMatLLL(mat);
  m1 := MutableTransposedMat(HermiteIntMatLLL(MutableTransposedMat(m1)));
  if not IsDiagonalMat(m1) then
    return SmithIntMatLLL(m1);
  fi;
  d :=  List([1..Minimum(Length(m1), Length(m1[1]))], i -> m1[i][i]);
  for i in [1..Length(d)] do
    for j in [i+1..Length(d)] do
      g := GcdInt(d[i],d[j]);
      if g<>d[i] then
        d[j] := d[i]*d[j]/g;
        d[i] := g;
      fi;
    od;
  od;
  
  for i in [1..Length(d)] do
    m1[i][i] := d[i];
  od;
  
  return m1;
end);

InstallGlobalFunction(SmithIntMatLLLTrans, function(mat)
  local   m1,  m2,  L,  R,  rowL,  colR,  c,  k,  d,  tmp,  g, x;
  # first Hermite normal form
  m1 := HermiteIntMatLLLTrans(mat);
  # Hermite normal form of the transposed of the Hermite normal form
  m2 := HermiteIntMatLLLTrans(MutableTransposedMat(m1[1]));
  # start again if not yet diagonal
  if not IsDiagonalMat(m2[1]) then
    x := SmithIntMatLLLTrans(MutableTransposedMat(m2[1]));
    return [x[1], x[2]*m1[2], MutableTransposedMat(m2[2])*x[3]];
  fi;
  # make diagonal entries dividing the next ones
  L := m1[2];
  R := m2[2];
  # a terrible "Immutable" is necessary here
  mat := MutableTransposedMat(m2[1]);
  c := [  ];
  for k  in [1..Minimum(Length(mat), Length(mat[1]))] do
    if mat[k][k] <> 0  then
      Add( c, mat[k][k] );
    fi;
  od;
  for d  in [1..Length( c )]  do
    for k  in [d+1..Length(c)] do
      if c[k] mod c[d] <> 0  then
        tmp := GcdRepresentation( c[d], c[k] );
        g := tmp * [ c[d], c[k] ];
        rowL := L[d];
        L[d] := rowL + tmp[2] * L[k];
        L[k] := -((1 - c[k] / g * tmp[2]) * L[k] - c[k] / g * rowL);
        colR := R[d];
        R[d] := R[k] + tmp[1] * colR;
        R[k] := (1 - c[d] / g * tmp[1]) * colR - c[d] / g * R[k];
        c[k] := c[k] * c[d] / g;
        c[d] := g;
      fi;
    od;
    mat[d][d] := c[d];
  od;

  return [mat, L, MutableTransposedMat(R)]; 
end);
