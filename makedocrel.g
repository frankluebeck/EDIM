##  this creates the documentation, needs: GAPDoc package, latex, pdflatex,
##  mkindex, dvips
##  

RequirePackage("GAPDoc");

MakeGAPDocDoc("doc", "edim", [], "EDIM", "../../..", "MathJax");
# need another run
Exec("cd doc; pdflatex edim.tex >> /dev/null; mv edim.pdf manual.pdf");

GAPDocManualLab("EDIM");

CopyHTMLStyleFiles("doc");

QUIT;

