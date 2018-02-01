#############################################################################
##  
##  PackageInfo.g file for the   EDIM   package.               Frank Lübeck
##  

SetPackageInfo( rec(

##  This is case sensitive, use your preferred spelling.
PackageName := "EDIM",

##  See '?Extending: Version Numbers' in GAP help for an explanation
##  of valid version numbers.
Version := "1.3.3",

##  Release date of the current version in dd/mm/yyyy format.
Date := "30/01/2018",

##  Subtitle, a header line for the package, for example a long form of the 
##  package name. (Should fit on one line together with the package name.)
Subtitle := "Elementary Divisors of Integer Matrices",

##  URL of the archive(s) of the current package release, but *without*
##  the format extension(s), like '.tar.gz' or '-win.zip', which are given next.
##  The archive file name must be changed with each version of the archive
##  (and probably somehow contain the package name an version).
ArchiveURL := "http://www.math.rwth-aachen.de/~Frank.Luebeck/EDIM/edim-1.3.3",

##  All provided formats as list of file extensions, separated by white
##  space or commas.
##  Currently recognized formats are:
##      .tar.gz    the UNIX standard
##      .tar.bz2   compressed with 'bzip2', often smaller than with gzip
##      -win.zip   zip-format for DOS/Windows, text files must have DOS 
##                 style line breaks (CRLF)
##  
##  In the future we may also provide .deb or .rpm formats which allow
##  a convenient installation and upgrading on Linux systems.
##  
ArchiveFormats := ".tar.bz2 .tar.gz -win.zip",


##  Information about authors and maintainers. Specify for each person a 
##  record with the following information:
##  
##     rec(
##     # these are compulsory, characters are interpreted as latin-1, so
##     # German umlauts and other western European special characters are ok:
##     LastName := "Müller",
##     FirstNames := "Fritz Eduard",
##  
##     # At least one of the following two entries must be given and set 
##     # to 'true' (an entry can be left out if value is not 'true'):
##     IsAuthor := true;
##     IsMaintainer := true;
##  
##     # At least one of the following three entries must be given.
##     # - preferably email address and WWW homepage
##     # - postal address not needed if email or WWW address available
##     # - if no contact known, specify postal address as "no address known"
##     Email := "Mueller@no.org",
##     # complete URL, starting with protocol
##     WWWHome := "http://www.no.org/~Mueller",
##     # separate lines by '\n'
##     PostalAddress := "Dr. F. Müller\nNo Org Institute\nNo Place 13\n\
##     12345 Notown\nNocountry"
##     
##     # If you want, add one or both of the following entries:
##     Place := "Notown",
##     Institution := "Institute for Nothing"
##     )
##  
Persons := [
  rec(
  LastName := "Lübeck",
  FirstNames := "Frank",
  IsAuthor := true,
  IsMaintainer := true,
  Email := "Frank.Luebeck@Math.RWTH-Aachen.De",
  WWWHome := "http://www.math.rwth-aachen.de/~Frank.Luebeck",
  PostalAddress := "Dr. Frank Lübeck\nLehrstuhl D für Mathematik\nRWTH Aachen\nTemplergraben 64\n52062 Aachen\nGERMANY\n",
  Place := "Aachen",
  Institution := "Lehrstuhl D für Mathematik, RWTH Aachen"
  )
],

##  Status information. Currently the following cases are recognized:
##    "accepted"      for successfully refereed packages
##    "deposited"     for packages for which the GAP developers agreed 
##                    to distribute them with the core GAP system
##    "dev"           for development versions of packages 
##    "other"         for all other packages
##
Status := "accepted",

##  You must provide the next two entries in case of status "accepted":
# format: 'name (place)'
CommunicatedBy := "Mike Atkinson (St Andrews)",
# format: mm/yyyy
AcceptDate := "08/1999",

##  For a central overview of all packages and a collection of all package
##  archives it is necessary to have two files accessible which should be
##  contained in each package:
##     - A README file, containing a short abstract about the package
##       content and installation instructions.
##     - The file you are currently reading or editing!
##  You must specify URLs for these two files, these allow to automate 
##  the updating of package information on the GAP Website, and inclusion
##  and updating of the package in the GAP distribution.
##  
README_URL := "http://www.math.rwth-aachen.de/~Frank.Luebeck/EDIM/README",
PackageInfoURL := 
           "http://www.math.rwth-aachen.de/~Frank.Luebeck/EDIM/PackageInfo.g",

##  Here you  must provide a short abstract explaining the package content 
##  in HTML format (used on the package overview Web page) and an URL 
##  for a Webpage with more detailed information about the package
##  (not more than a few lines, less is ok):
# Please, use '<span class="pkgname">GAP</span>' and
# '<span class="pkgname">MyPKG</span>' for specifing package names.
AbstractHTML := "This package provides  a collection of functions for \
computing the Smith normal form of integer matrices and some related \
utilities.",

PackageWWWHome := "http://www.math.rwth-aachen.de/~Frank.Luebeck/EDIM",

##  On the GAP Website there is an online version of all manuals in the
##  GAP distribution. To handle the documentation of a package it is
##  necessary to have:
##     - an archive containing the package documentation (if possible
##       HTML and PDF-format)
##     - the start file of the HTML documentation (if provided), relative to
##       package root
##     - the PDF-file relative to the package root
##  For links to other package manuals or the GAP manuals one can assume 
##  relative paths as in a standard GAP installation. 
##  Also, provide the information which is currently given in your packages 
##  init.g file in the command DeclarePackage(Auto)Documentation
##  (for future simplification of the package loading mechanism).
##  
##  Please, remove all unnecessary files (.log, .aux, .dvi, .ps, ...) from
##  the documentation archive.
##  
# in case of several help books give a list of entries here:
PackageDoc := rec(
  # use same as in GAP            
  BookName := "EDIM",
  # a list of directories and/or file names in the package archive which
  # contains the documentation
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0.html",
  PDFFile := "doc/manual.pdf",
  # the path to the .six file used by GAP's help system
  SixFile := "doc/manual.six",
  # a longer title of the book, this together with the book name should
  # fit on a single text line (appears with the '?books' command in GAP)
  LongTitle := "Elementary Divisors of Integer Matrices",
  # Should this help book be autoloaded when GAP starts up? This should
  # usually be 'true', otherwise say 'false'. 
  Autoload := true
),


##  Are there restrictions on the operating system for this package? Or does
##  the package need other packages to be available?
Dependencies := rec(
  # GAP version, use version strings for specifying exact versions,
  # prepend a '>=' for specifying a least version.
  GAP := "4.5",
  # list of pairs [package name, (least) version],  package name is case
  # insensitive, least version denoted with '>=' prepended to version string.
  # without these, the package will not load
  NeededOtherPackages := [["GAPDoc", ">= 1.5"]],
  # without these the package will issue a warning while loading
  SuggestedOtherPackages := [],
  # needed external conditions (programs, operating system, ...)  provide 
  # just strings as text or
  # pairs [text, URL] where URL  provides further information
  # about that point.
  # (no automatic test will be done for this, do this in your 
  # 'AvailabilityTest' function below)
  ExternalConditions := ["UNIX for the kernel function 'ElementaryDivisorsPPartRkExpSmall'"]
                      
),

#  Provide a test function for the availability of this package, see
#  documentation of 'Declare(Auto)Package', this is the <tester> function.
#  For packages which will not fully work, use 'Info(InfoWarning, 1,
#  ".....")' statements. For packages containing nothing but GAP code,
#  just say 'ReturnTrue' here.
#  (When this is used for package loading in the future the availability
#  tests of other packages, as given above, will be done automatically and
#  need not be included here.)
AvailabilityTest := function()
  if not "ediv" in SHOW_STAT() and 
     Filename(DirectoriesPackagePrograms("edim"), "ediv.so") = fail then
    LogPackageLoadingMessage( PACKAGE_WARNING,
              [ "The EDIM kernel function 'ElementaryDivisorsPPartRkExpSmall'",
                "is not available.",
                "It is recommended to compile this function, see",
                "'?Installation of the EDIM package'" ] );
  fi;
  return true;
end,

Autoload := false,

# we now use the default banner
##  BannerString := Concatenation(
##  "    ######################################################################\n",
##  "    ##                                                                  ##\n",
##  "    ##       EDIM ",
##  ~.Version, " (Elementary Divisors and Integer Matrices)      ##\n",
##  "    ##                                                                  ##\n",
##  "    ##   Questions and remarks to: Frank.Luebeck@Math.RWTH-Aachen.De    ##\n",
##  "    ##                                                                  ##\n",
##  "    ######################################################################\n\n"
##  ),

##  Optional, but recommended: path relative to package root to a file which 
##  contains as many tests of the package functionality as sensible.
#TestFile := "tst/testall.g",

##  Optional: Here you can list some keyword related to the topic 
##  of the package.
Keywords := ["Smith normal form", "p-adic", "rational matrix inversion"]

));


