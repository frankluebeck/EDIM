/***************************************************************************
**
*A  ediv.c               EDIM-mini-package                     Frank Lübeck
**
**  
*Y  Copyright (C) 1999  Lehrstuhl D f\"ur Mathematik, RWTH Aachen
**  
**  `ElementaryDivisorsPPartRkExpSmall' as kernel function.
**  
*/


/* read GAP source header files with a combined header file */

#include        "src/compiled.h"          /* GAP headers                */


#include        <stdio.h>
#include        <limits.h>

/*          the functions         */

/* a helper function for computing k^-1 mod p */
long invmodpcint(long k, long p)
{
  long f, g, h, fk, gk, hk, q;
  if (0 <= k) {
    f = k;
    fk = 1;
  }
  else {
    f = -k;
    fk = -1;
  }
  g = p;
  gk = 0;
  while (g != 0) {
    q =  f/g;
    h = g;
    hk = gk;
    g = f - q * g;
    gk = fk - q * gk;
    f = h;
    fk = hk;
  }
  return fk % p;
}

/* This function is transformed from a stand alone  C-program, which looks
   very much like the corresponding GAP-function.
   Instead of 'malloc' we use here 'NewBag' with type 'T_DATOBJ'. */ 

Obj FuncElementaryDivisorsPPartRkExpSmall(
					 Obj self,
					 Obj A,
					 Obj pobj,
					 Obj rkobj,
					 Obj robj,
					 Obj ilobj) /* info level */
{
  unsigned long m, n, nn, r, rk, ch, chmax, p, pr, i, j, 
                ii, i0, i1, x, c, pos, il;
  unsigned long *A1, *A2, *Tp, *res, *inv, *vv;
  Obj A1obj, A2obj, Tpobj, resobj, invobj, probj, obj, row;

  /*  change args  */
  if (! (IS_PLIST(A) || IS_POSOBJ(A)) || LEN_PLIST(A) < 1 || 
      ! (IS_PLIST(A) || IS_POSOBJ(ELM_PLIST(A, 1)))) {
     ErrorQuit("A must be an integer matrix",0L,0L); 
  }
  m = (unsigned long) LEN_PLIST(A);
  n = (unsigned long) LEN_PLIST(ELM_PLIST(A, 1));
  if (! IS_INTOBJ(pobj)) { 
     ErrorQuit("p must be a small integer (not a %s)",(Int)TNAM_OBJ(pobj),0L); 
  }
  p = (unsigned long) INT_INTOBJ(pobj);
  if (! IS_INTOBJ(rkobj)) { 
     ErrorQuit("rk must be a small integer (not a %s)",(Int)TNAM_OBJ(rkobj),0L);
  }
  rk = (unsigned long) INT_INTOBJ(rkobj);
  if (! IS_INTOBJ(robj)) { 
     ErrorQuit("r must be a small integer (not a %s)",(Int)TNAM_OBJ(robj),0L);
  }
  r = (unsigned long) INT_INTOBJ(robj);
  if (! IS_INTOBJ(ilobj)) { 
     ErrorQuit("il must be a small integer (not a %s)",(Int)TNAM_OBJ(ilobj),0L);
  }
  il = (unsigned long) INT_INTOBJ(ilobj);

  /* pr = p^(r+1) */
  probj = PowInt(pobj, SumInt(robj, INTOBJ_INT(1)));
  if (! IS_INTOBJ(probj)) { 
     ErrorQuit("exponent too large, see ?ElementaryDivisorsPPartRkExpSmall",0L,0L);
  }
  /* p^(r+1)-1 must fit at least (p-1) times into an unsigned long */
  pr = (unsigned long) INT_INTOBJ(probj);
  if (ULONG_MAX/(pr-1) < (p-1)) {
     ErrorQuit("exponent too large, see ?ElementaryDivisorsPPartRkExpSmall",0L,0L);
  }
  /* max sum of coeffs of numbers of size (p^(r+1)-1) before reduction is 
     necessary to avoid integer overflow */
  chmax = ULONG_MAX/(pr-1)-1;
  A2obj = NewBag(T_DATOBJ, (n+1)*(m+1)*sizeof(unsigned long));
  A2 = (unsigned long *)ADDR_OBJ(A2obj);
  A2[0] = m;
  nn = n+1;
  /* reduce matrix entries modulo p^(r+1) */
  for (i=1; i<=m; i++) {
    row = ELM_PLIST(A, i);
    if (! IS_PLIST(row) || LEN_PLIST(row) < n) {
       ErrorQuit("A must be an integer matrix",0L,0L);
    }
    for (j=1; j<=n; j++) {
      obj = ELM_PLIST(ELM_PLIST(A, i), j);
      if (!(IS_INTOBJ(obj) || TNUM_OBJ(obj)==T_INTPOS || TNUM_OBJ(obj)==T_INTNEG)) {
         ErrorQuit("matrix entry must be integer (not a %s)",
                                              (Int)TNAM_OBJ(obj),0L);
      }
      if (LtInt(probj, obj) || LtInt(obj, INTOBJ_INT(0))) {
	obj = ModInt(obj, probj);
        /* this could have changed because of a garbage collection */
	A2 = (unsigned long *)ADDR_OBJ(A2obj);
	A2[i*nn+j] = (unsigned long) INT_INTOBJ(obj);
      }
      else {
	A2 = (unsigned long *)ADDR_OBJ(A2obj);
	A2[i*nn+j] = (unsigned long) INT_INTOBJ(obj);
      }
    }
  }

  /* allocating space for local variables */
  A1obj = NewBag(T_DATOBJ, (n+1)*(m+1)*sizeof(unsigned long));
  invobj = NewBag(T_DATOBJ, (n+1)*sizeof(unsigned long));
  resobj = NewBag(T_DATOBJ, (r+2)*sizeof(unsigned long));
  Tpobj = NewBag(T_DATOBJ, (n+1)*sizeof(unsigned long));
  A1 = (unsigned long *)ADDR_OBJ(A1obj);
  A1[0] = 0UL;
  res = (unsigned long *)ADDR_OBJ(resobj);
  res[0] = 0UL;
  inv = (unsigned long *)ADDR_OBJ(invobj);
  Tp = (unsigned long *)ADDR_OBJ(Tpobj);
  A2 = (unsigned long *)ADDR_OBJ(A2obj);

  /* from now on the pointers above are safe: we only  manipulate  the data
     in the allocated bags and don't use any GAP function which could cause
     a garbage collection */
  for (j=1; j<=n; Tp[j]=j, j++);

  while (A1[0]<rk && r+1>0) {
    if (il>0) {fprintf(stderr, "#"); fflush(stderr);}
    i0 = A2[0];
    A2[0] = 0;
    for (ii=1; ii<=i0; ii++) {
      vv = A2+ii*nn;
      /* ch counts how many times p^(r+1) the biggest possible entry of
         v is */
      for (i=1, ch=1; i<=A1[0]; i++) {
	c = ((vv[Tp[i]] % p)*(p-inv[i])) % p;
	if (c != 0) {
          ch += c;
	  if (ch>=chmax){
	    /* reduce the entries */
	    for (j=1; j<=n; vv[j] %= pr, j++);
	    ch = c+1;
	  }
	  for (j=1; j<=n; j++) {
	    vv[j] += (c * A1[i*nn+j]);
	  }
	}
      }
      for (j=1; j<=n; vv[j] %= pr, j++);
      for (j=pos=A1[0]+1; (pos<=n) && (vv[Tp[pos]]%p == 0); pos++);
      if (pos>n) {
	if (il>0) {fprintf(stderr, "-");fflush(stderr);}
	A2[0]++;
	i1 = A2[0];
	for (j=1; j<=n; j++) {
	  A2[i1*nn+j] = vv[j]/p;
	}
      }
      else {
	if (il>0) {fprintf(stderr, "+");fflush(stderr);}
	A1[0]++;
	i1 = A1[0];
	for (j=1; j<=n; j++) {
	  A1[i1*nn+j] = vv[j];
	}
	if (pos != i1) {
	    x = Tp[i1];
	    Tp[i1] = Tp[pos];
	    Tp[pos] = x;
	}
	inv[i1] = invmodpcint(A1[i1*nn+Tp[i1]] % p, p);
      }
    }
    res[0]++;
    if (il>0) {fprintf(stderr, "\n#Rank found: %ld\n", A1[0]); fflush(stderr);}
    res[res[0]] = rk-A1[0];
    r--;
  }
  if (A1[0]==rk) { 
    RetypeBag(resobj, T_PLIST);
    i0 = res[0];
    SET_LEN_PLIST(resobj, i0);
    for (i=1; i<=i0; i++) {
      SET_ELM_PLIST(resobj, i, INTOBJ_INT(res[i]));
    }
  }
  else {
    if (il>0) {
      fprintf(stderr, "#exponent too small or rank too big. \n"); 
      fflush(stderr);
    }
    resobj = Fail;
  }
  return resobj;
}

/*F * * * * * * * * * * * * * initialize package * * * * * * * * * * * * * * *
*/



/****************************************************************************
**

*V  GVarFilts . . . . . . . . . . . . . . . . . . . list of filters to export
*/
/*static StructGVarFilt GVarFilts [] = {

    { "IS_BOOL", "obj", &IsBoolFilt,
      IsBoolHandler, "src/bool.c:IS_BOOL" },

    { 0 }

};   ?????*/

/****************************************************************************
**

*V  GVarFuncs . . . . . . . . . . . . . . . . . . list of functions to export
*/
static StructGVarFunc GVarFuncs [] = {

  { "ElementaryDivisorsPPartRkExpSmall", 5, "A, p, rk, r, il", 
    FuncElementaryDivisorsPPartRkExpSmall, 
    "ediv.c:ElementaryDivisorsPPartRkExpSmall" },

  { 0 }

};



/**************************************************************************

*F  InitKernel( <module> )  . . . . . . . . initialise kernel data structures
*/
static Int InitKernel (
    StructInitInfo *    module )
{

    /* init filters and functions                                          */
    InitHdlrFuncsFromTable( GVarFuncs );

    /* return success                                                      */
    return 0;
}


/****************************************************************************
**
*F  InitLibrary( <module> ) . . . . . . .  initialise library data structures
*/
static Int InitLibrary (
    StructInitInfo *    module )
{
  /*    UInt            gvar;
	Obj             tmp; */

    /* init filters and functions                                          */
    /* printf("Init El..Small\n");fflush(stdout); */
    InitGVarFuncsFromTable( GVarFuncs );

    /* return success                                                      */
    return 0;
}


/****************************************************************************
**
*F  InitInfopl()  . . . . . . . . . . . . . . . . . table of init functions
*/
/* <name> returns the description of this module */
static StructInitInfo module = {
#ifdef EDIVSTATIC
    .type = MODULE_STATIC,
#else
    .type = MODULE_DYNAMIC,
#endif
    .name = "ediv",
    .initKernel = InitKernel,
    .initLibrary = InitLibrary,
};

#ifndef EDIVSTATIC
StructInitInfo * Init__Dynamic ( void )
{
 return &module;
}
#endif

StructInitInfo * Init__ediv ( void )
{
  return &module;
}

