LoadPackage( "EDIM" );

TestDirectory( DirectoriesPackageLibrary("EDIM", "tst"),
            rec(exitGAP := true ) );

# Should never get here
FORCE_QUIT_GAP(1);
