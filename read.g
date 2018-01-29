#############################################################################
##
#A  read.g               EDIM-mini-package                     Frank Lübeck
##
##  
#Y  Copyright (C) 1999  Lehrstuhl D f\"ur Mathematik, RWTH Aachen
##  
##  Reading the library of the package.
##  

ReadPackage("edim", "lib/util.gi");
ReadPackage("edim", "lib/edim.gi");
ReadPackage("edim", "lib/hmmlll.gi");

# redefine the Info handler and output, we use the plain one from GAPDoc
SetInfoHandler(InfoEDIM, PlainInfoHandler);
SetInfoOutput(InfoEDIM, "*errout*");

