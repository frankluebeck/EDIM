
LoadPackage("EDIM");

exmpls := ExtractExamples(".", "edim.xml", [], "all");

RunExamples(exmpls);
#RunExamples(exmpls, rec(changeSources := true));

