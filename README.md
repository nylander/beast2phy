beast2phy -- Convert output trees from BEAST to other formats  
=============================================================

Convert output trees from BEAST to other formats

What:

    beast2phy.pl with helper scripts

Version:

    04/15/2008 01:53:14 PM CEST
    04/26/2011 12:50:44 PM CEST

By:

    Johan Nylander

Usage:

    beast2phy.pl [--format=<nex|altnex|mb|phy> --decimals=N --outfile=OUTFILE]  INFILE

    beast2mb INFILE

    beast2altnex INFILE

    beast2newick INFILE


Description:

    Scripts to convert the output from BEAST/TreeAnnotator to formats readable by other software.


Files:

    beast2phy.pl    -- the main Perl script for MCMC sampled trees. Try 'beast2phy.pl --help'

    beast2mb        -- beast to "strict" MrBayes output

    beast2altnex    -- beast to alt nexus (no translation table)

    beast2newick    -- beast to newick (phylip) format

    trees.txt       -- example beast file

    trees2.txt      -- example beast file

    trees.con       -- example beast file


Installation:

    Put all files in the PATH (e.g., /usr/local/bin)


Warning:

    The script can be VERY slow on large files (e.g., minutes on a 250 MB tree file)

