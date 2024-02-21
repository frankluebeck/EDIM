#############################################################################
##  
##  PackageInfo.g file for the   EDIM   package.               Frank Lübeck
##  
SetPackageInfo( rec(

PackageName := "EDIM",
Version := "1.3.8",
##  dd/mm/yyyy 
Date := "21/02/2024",
License := "GPL-2.0-or-later",
Subtitle := "Elementary Divisors of Integer Matrices",
# without extension
ArchiveURL := "https://www.math.rwth-aachen.de/~Frank.Luebeck/EDIM/EDIM-1.3.8",
ArchiveFormats := ".tar.bz2  .tar.gz   -win.zip",
SourceRepository := rec(Type := "git", 
                        URL := "https://github.com/frankluebeck/EDIM" ),
Persons := [
  rec(
  LastName := "Lübeck",
  FirstNames := "Frank",
  IsAuthor := true,
  IsMaintainer := true,
  Email := "Frank.Luebeck@Math.RWTH-Aachen.De",
  WWWHome := "https://www.math.rwth-aachen.de/~Frank.Luebeck",
  PostalAddress := "Dr. Frank Lübeck\nLehrstuhl D für Mathematik\nRWTH Aachen\nTemplergraben 64\n52062 Aachen\nGERMANY\n",
  Place := "Aachen",
  Institution := "Lehrstuhl D für Mathematik, RWTH Aachen"
  )
],
Status := "accepted",
CommunicatedBy := "Mike Atkinson (St Andrews)",
# mm/yyyy
AcceptDate := "08/1999",
README_URL := "https://www.math.rwth-aachen.de/~Frank.Luebeck/EDIM/README",
PackageInfoURL := 
           "https://www.math.rwth-aachen.de/~Frank.Luebeck/EDIM/PackageInfo.g",
AbstractHTML := "This package provides  a collection of functions for \
computing the Smith normal form of integer matrices and some related \
utilities.",
PackageWWWHome := "https://www.math.rwth-aachen.de/~Frank.Luebeck/EDIM",
PackageDoc := rec(
  BookName := "EDIM",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0.html",
  PDFFile := "doc/manual.pdf",
  SixFile := "doc/manual.six",
  LongTitle := "Elementary Divisors of Integer Matrices",
  Autoload := true
),
Dependencies := rec(
  GAP := "4.12",
  NeededOtherPackages := [["GAPDoc", ">= 1.6"]],
  SuggestedOtherPackages := [],
  ExternalConditions := 
        ["UNIX for the kernel function 'ElementaryDivisorsPPartRkExpSmall'"]
),

AvailabilityTest := function()
  if not IsKernelExtensionAvailable("EDIM","ediv") then
    LogPackageLoadingMessage( PACKAGE_WARNING,
              [ "The EDIM kernel function 'ElementaryDivisorsPPartRkExpSmall'",
                "is not available.",
                "It is recommended to compile this function, see",
                "'?Installation of the EDIM package'" ] );
  fi;
  return true;
end,
Autoload := false,
TestFile := "tst/edim.tst",
Keywords := ["Smith normal form", "p-adic", "rational matrix inversion"]
));


