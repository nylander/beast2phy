#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  beast2phy.pl
#
#        USAGE:  ./beast2phy.pl --format=nex|altnex|mb|phy --decimals=N --outfile=OUTFILE  INFILE
#
#  DESCRIPTION:  Tries to convert a Beast .con file (output from TreeAnnotator)
#                to a tree in Phylip (Newick) format (readable by, e.g. read.tree()
#                in APE).
#                Does also convert the file with MCMC sampled trees from BEAST.
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#    BUGS/TODO:  branch lengths written in scientific notation?
#        NOTES:  ---
#       AUTHOR:  Johan A. A. Nylander (JN), <jnylander @ users.sourceforge.net>
#      COMPANY:  SU
#      VERSION:  1.0
#      CREATED:  09/13/2008 06:17:14 PM CEST
#     REVISION:  04/26/2011 12:35:55 PM CEST
#===============================================================================

use strict;
use warnings;
use Getopt::Long;



#===  FUNCTION  ================================================================
#         NAME:  find_taxa_block
#      VERSION:  04/15/2009 09:48:19 AM CEST
#  DESCRIPTION:  finds start and end of taxa block
#   PARAMETERS:  ???
#      RETURNS:  "0"/"1", $start, $end, 
#         TODO:  very ad hoc solutions overall. rewrite?
#===============================================================================
sub find_taxa_block {

    my ($infile) = (@_);

    my $block_start = q{};
    my $block_end   = q{};
    my $has_block   = q{};

    open my $INFILE, "<", $infile or die "in find_taxa_block: could not open infile; $! \n";

    while(<$INFILE>) {
        if (/^\s*End\s*;/i) {
            $block_end = $.;
            last;
        }
        elsif (/^\s*Begin\s+taxa\s*;/i) {
            $block_start = $.;
        }
    }

    close($INFILE) or warn "in find_taxa_block: could not close file handle: $! \n";

    if ($block_start eq '') {
        $has_block = 0;
        $block_start = '';
        $block_end = '';
    }
    else {
        $has_block = 1;
    }

    return($has_block, $block_start, $block_end);

} # end of find_taxa_block


#===  FUNCTION  ================================================================
#         NAME:  remove_comments
#      VERSION:  09/13/2008 05:21:32 PM CEST
#  DESCRIPTION:  removes comments (text withing square brackets) from a tree description
#   PARAMETERS:  tree string
#      RETURNS:  tree string
#         TODO:  Stream line this function and or the call to/from this function! 
#                This is where the script spend most of it's time
#===============================================================================
sub remove_comments {

    my ($tree) = @_;

    my $outside = 1;
    my $inside = 0;
    my (@array) = split //, $tree;
    my (@return) = ();
    my $tree_string = q{};

    foreach my $char (@array) {
        if($char eq '[') {
            $inside = 1;
            $outside = 0;
            next;
        }
        elsif ($char eq ']') {
            $outside = 1;
            $inside = 0;
            next;
        }
        elsif ($outside) {
            push @return, $char;
        }
    }
    
    $tree_string = join '', @return;

    return($tree_string);
 
} # end of remove_comments


#===  FUNCTION  ================================================================
#         NAME:  round_numbers
#      VERSION:  09/13/2008 05:20:49 PM CEST
#  DESCRIPTION:  Round numbers
#   PARAMETERS:  ndecimals, tree string
#      RETURNS:  tree string
#         TODO:  Handle cases where printing with too few decimals gives 0-length branches
#===============================================================================
sub round_numbers {

    my ($ndecimals, $tree) = @_;
    
    $tree =~ s/(\d+\.[\d|E-]+)/sprintf("%.${ndecimals}f", $1)/gei;

    return($tree);
 
} # end of round_numbers


#===  FUNCTION  ================================================================
#         NAME:  read_labels
#      VERSION:  04/26/2011 12:34:11 PM CEST
#  DESCRIPTION:  reads sequence (taxon) labels and associates them with numbers
#                The function tries to read within the trees block.
#   PARAMETERS:  file name
#      RETURNS:  Hash with taxon translation table
#         TODO:  ??? 
#===============================================================================
sub read_labels {

    my ($file) = @_;
    my %hash = ();
    my $DEBUG = 0;

    open (my $FILE, '<', $file) or die "$0 : failed to open input file $file : $!\n";
    while(<$FILE>) {
        my $line = $_;
        chomp($line);
        if($line =~ m/^\s*tree\s+/i) {
            last;
        }
        elsif($line =~ m/^\s+(\d+)\s+([\w|\s|'|\d|,|;]+)\Z/) { # capture the number, and the
            my $number = $1;                                   # taxon name allowing for single-
            my $name = $2;                                     # quoted taxon names
            $name =~ s/,\Z//;
            $name =~ s/;\Z//;
            $hash{$number} = $name;
        }
    }
    close ($FILE) or warn "$0 : failed to close file $file : $!\n";

    if ($DEBUG) {
        print STDERR "\nhash in read_labels:\n";
        for my $key ( keys %hash ) {
            my $value = $hash{$key};
            print STDERR "$key => $value\n";
        }
    }

    return %hash;

} # end of read_labels


#===  FUNCTION  ================================================================
#         NAME:  replace_taxon_numbers
#      VERSION:  02/02/2007 01:48:44 AM CET
#  DESCRIPTION:  replaces numbers with sequence (taxon) labels
#   PARAMETERS:  tree string and reference to hash holding translation table.
#                Uses global variable $infile
#      RETURNS:  tree string with seq labels instead of numbers
#         TODO:  ?
#===============================================================================
sub replace_taxon_numbers {

    my ($tree, $hash_reference) = @_;

    foreach my $number (keys %$hash_reference) { # $key is number, $value is label
        my $label = $hash_reference->{$number};
        if ($tree =~ /,$number:/) {
            $tree =~ s/,$number:/,$label:/;      # ',123:' => ',foo:'
        }
        elsif ($tree =~ /\($number:/) {
            $tree =~ s/\($number:/\($label:/;    # '(123:' => '(foo:'
        }
        elsif ($tree =~ /\($number,/) {
            $tree =~ s/\($number,/\($label,/;    # '(123,' => '(foo,'
        }
        elsif ($tree =~ /,$number,/) {
            $tree =~ s/,$number,/,$label,/;      # ',123,' => ',foo,'
        }
        else {
            $tree =~ s/,$number\)/,$label\)/;    # ',123)' => ',foo)'
        }
    }

    return $tree;

} # end of replace_taxon_numbers


#===  FUNCTION  ================================================================
#         NAME:  replace_tree_name
#      VERSION:  01/21/2009 12:13:51 PM CET
#  DESCRIPTION:  replaces tree names with names readable by mrbayes
#   PARAMETERS:  
#      RETURNS:  tree string with replaced tree name
#         TODO:  ?
#===============================================================================
sub replace_tree_name {

    my ($line)      = @_;
    my @tree_pieces = ();
    my $tree        = q{};
    my $space       = q{ };

    my @pieces = split /\s+/, $line;

    foreach my $piece (@pieces) {
        $piece =~ s/^STATE_(\d+)$/rep\.$1/;
        push @tree_pieces, $piece; 
    }
    $tree = join $space, @tree_pieces; 
    
    $tree = "$tree\n"; # hack

    return $tree;

} # end of replace_tree_name


#===  FUNCTION  ================================================================
#         NAME:  remove_tree_name
#      VERSION:  09/13/2008 03:04:10 PM CEST
#  DESCRIPTION:  Removes the tree name
#   PARAMETERS:  Line with tree string. Uses global variable $infile
#      RETURNS:  Tree string
#         TODO:  The regexp will not work correctly if there is a space between
#                the last closing bracket and the trailing ';'. Could be solved
#                like this:
#                if ( $piece =~ ^\(.+\)$ ) { # starts with ''( and ends in ')'
#                    $tree = $piece . ";";
#                }
#                elsif ($piece =~ ^\(.+[\);]$) { # starts with '(' and ends in ');'
#                    $tree = $piece . ";";
#                }
#                Any spaces in the tree description will of course create
#                truncated trees, however.
#===============================================================================
sub remove_tree_name {
    
    my ($line) = @_;
    my $tree   = '';

    my @pieces = split /\s+/, $line;

    foreach my $piece (@pieces) {
        if ($piece =~ /^\(.+[;]$/) { # if piece starts with '(' and ends in ';'
            $tree = "$piece\n";        # hack: needed to add line break
        }
    }

    if ($tree eq '') {
        die "Warning: Could not read tree format.\nAcceptable format: tree name = (1,(2,3));\n";
    }

    return $tree;

} # end of remove_tree_name


#===  FUNCTION  ================================================================
#         NAME:  MAIN 
#      VERSION:  04/15/2009 10:59:07 AM CEST
#  DESCRIPTION:  ???
#   PARAMETERS:  ???
#      RETURNS:  ???
#         TODO:  ???
#===============================================================================
MAIN:

## Globals
my $DEBUG             = 0;   # 0 or 1
my $outfile           = q{}; # Empty
my $infile            = q{};
my $format            = q{};
my $nexus             = q{};
my $altnexus          = q{};
my $mb                = q{};
my $phy               = q{};
my $using_outfile     = q{};
my $has_taxa_block    = q{};
my $taxa_block_start  = q{};
my $taxa_block_end    = q{};
my $space             = q{};
my $decimals          = q{}; # Number of decimals if PPs and branch lengths are to be rounded. Same nr applies to both.
my $decimals_default  = q{15};
my $PRINT_FH;
my %translation_table = ();
my $USAGE             = "--format=nex|altnex|mb|phy --decimals=N --outfile=OUTFILE INFILE\n\n";


## Arguments
if (@ARGV < 1) {
    print "\n Usage:\n\n $0 $USAGE";
    exit(0);
}
else {
    GetOptions('help'       => sub { print "\n Usage:\n\n $0 $USAGE"; exit(0); },
               'outfile:s'  => \$outfile,
               'decimals:i' => \$decimals,
               'format:s'   => \$format
              );
}

## Get infile
if ($infile eq '') {
    if (@ARGV > 0) { 
        $infile = shift(@ARGV);
        die "$0 : file $infile doesn't seem to exist: $!\n" unless (-e $infile);
    }
}
else {
    open (my $IN, '<', $infile) or die "$0 : failed to open input file $infile : $!\n";
}

## Set filehandle for printing.
if ($outfile eq '') {
    $PRINT_FH = *STDOUT; # Using the typeglob notation in order to use STDOUT as a variable
    $using_outfile = 0;
}
else {
    open ($PRINT_FH, '>', $outfile) or die "$0 : Failed to open output file $outfile : $!\n\n";
    $using_outfile = 1;
}

## Set output format (nexus, altnexus, mb, or phylip)
if ($format =~ /^a/i) {
    $altnexus = 1;
}
elsif ($format =~ /^n/i) {
    $nexus = 1;
}
elsif ($format =~ /^m/i) {
    $mb = 1;
}
elsif ($format =~ /^m/i) {
    $phy = 1;
}
else { # Default
    $phy = 1;
}

## Print file header
if($altnexus) {
    print $PRINT_FH "#NEXUS\n[Trees from file $infile]\nBegin trees;\n";
    $space = q{    }; # Four blanks
}
if($nexus) {
    print $PRINT_FH "#NEXUS\n[Trees from file $infile]\n";
    $space = q{    }; # Four blanks
}
if($mb) {
    print $PRINT_FH "#NEXUS\n[ID: 0123456789]\n";
    $space = q{   }; # 3 blanks
}

## Set number of decimals. Defaults to $decimals_default.
$decimals = $decimals_default unless ($decimals);

## Does infile have a taxa block?
($has_taxa_block, $taxa_block_start, $taxa_block_end) = find_taxa_block($infile);

## Open infile
open (my $IN, '<', $infile) or die "$0 : failed to open input file $infile : $!\n";

## Get the translation table
if ($altnexus or $phy) {
    %translation_table = read_labels($infile);
}

## Find the tree(s) and convert
while(<$IN>) {
    my $infile_line = $_;
    if (/^tree\s+/i) { # found tree string

        my $the_tree = $infile_line;
            if($DEBUG) {print STDERR "$the_tree\n"; warn "\n1 Made it here, hit return to continue\n" and getc();}

        $the_tree =~ s/posterior=(\d+.\d+)/\]$1\[/g; # replace pp with ']pp['
            if($DEBUG) {print STDERR "$the_tree\n"; warn "\n2 Made it here, hit return to continue\n" and getc();}

        $the_tree = remove_comments($the_tree);
            if($DEBUG) {print STDERR "$the_tree\n"; warn "\n3 Made it here, hit return to continue\n" and getc();}

        $the_tree = round_numbers($decimals, $the_tree);
            if($DEBUG) {print STDERR "$the_tree\n"; warn "\n4 Made it here, hit return to continue\n" and getc();}

        if($mb) {
            $the_tree = replace_tree_name($the_tree);
            if($DEBUG) {print STDERR "$the_tree\n"; warn "\n5 Made it here, hit return to continue\n" and getc();}
        }

        if ($phy) {
            $the_tree = remove_tree_name($the_tree);
                if($DEBUG) {print STDERR "$the_tree\n"; warn "\n6 Made it here, hit return to continue\n" and getc();}
            $the_tree = replace_taxon_numbers($the_tree, \%translation_table);
                if($DEBUG) {print STDERR "$the_tree\n"; warn "\n7 Made it here, hit return to continue\n" and getc();}
        }

        if($altnexus) {
            $the_tree = replace_taxon_numbers($the_tree, \%translation_table);
                if($DEBUG) {print STDERR "$the_tree\n"; warn "\n8 Made it here, hit return to continue\n" and getc();}
        }
        print $PRINT_FH "$space", "$the_tree";
    }
    else {
        if ($phy or $altnexus) {
            # do nothing
        }
        else {
            if ($mb) {
                if ($has_taxa_block) {# handle taxa block!
                    if ($. < $taxa_block_start || $. > $taxa_block_end) {
                        print $PRINT_FH "$infile_line" unless ($infile_line =~ /^#|(^\s+$)/);
                    }
                }
                else {
                    print $PRINT_FH "$infile_line" unless ($infile_line =~ /^#|(^\s+$)/);
                }
            }
            else {
                if ($has_taxa_block) {# handle taxa block
                    if ($. < $taxa_block_start || $. > $taxa_block_end) {
                        print $PRINT_FH "$infile_line" unless ($infile_line =~ /^#/);
                    }
                }
                else {
                    print $PRINT_FH "$infile_line" unless ($infile_line =~ /^#/);
                }
            }
        }
    }
}

## Trailing end to altnexus
if ($altnexus) {
    print $PRINT_FH "END;\n\n";
}

## Close
close ($IN) or warn "$0 : failed to close input file $infile : $!\n";
if ($using_outfile) {
    close ($PRINT_FH) or warn "$0 : failed to close output file $outfile : $!\n";
}

__END__
