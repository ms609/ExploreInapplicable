use File::Find;

$ROOT = "C:/Research/ExploreInapplicable/";
$MATRIX_DIR = "matrices/";
$CHAR_TYPE_DIR = "charType/";
print " Finding matrices in ". $ROOT . $MATRIX_DIR;
find (\&char_template, $ROOT . $MATRIX_DIR);
print "\nDone.";

sub char_template() {
	if (-f and /^Wilso.*\.nex$/) {
    $nexus_filename    = $_;
    $raw_matrix_path   = $ROOT . $MATRIX_DIR    . $nexus_filename;
    $chartype_path     = $ROOT . $CHAR_TYPE_DIR . $nexus_filename;
    $chartype_path     =~ s/\.nex$/.txt/i;
    $extra_state_path  = $ROOT . 'extraState/'   . $nexus_filename;
    $ambiguous_path    = $ROOT . 'ambiguous/'    . $nexus_filename;
    $ambig_absent_path = $ROOT . 'ambigAbsent/'  . $nexus_filename;
    $inapplicable_path = $ROOT . 'inapplicable/' . $nexus_filename;
    
    print "\n- Processing $nexus_filename ... ";
    
    open (CHARTYPE, "<$chartype_path") or warn ("!! Can't open $chartype_path: $!\n");
    my @chartype = <CHARTYPE>;
    close CHARTYPE;
    foreach (@chartype) {
      if (length ($_) && !($_ =~ /^[INTX]/i)) {warn ("!! Invalid chartype!");};
    }
    
    open (MATRIX, "<$raw_matrix_path") or warn " ERROR: Can't open $raw_matrix_path.\n";
    @nexus = <MATRIX>;
    close MATRIX;
    
    open ($AMBIGUOUS , ">$ambiguous_path"   ) or warn "!! Can't open $ambiguous_path: $!\n";
    open ($AMBIGABS  , ">$ambig_absent_path") or warn "!! Can't open $ambig_absent_path: $!\n";
    open ($EXTRASTATE, ">$extra_state_path" ) or warn "!! Can't open $extra_state_path: $!\n";
    open ($INAPPLIC  , ">$inapplicable_path") or warn "!! Can't open $inapplicable_path: $!\n";
    my @matrix_files = ($AMBIGUOUS, $AMBIGABS, $EXTRASTATE, $INAPPLIC);
    
    # Write matrices for different analyses
    my @taxa = ();
    my $in_matrix = 0;
    foreach $line (@nexus) {
      $line_modified = 0;
      if ($line =~ /^\s*matrix\s*$/i) {
        $in_matrix = 1;
      } elsif ($line =~ /;/) {
        $in_matrix = 0;
      } elsif ($in_matrix && $line =~ /^(\s*[A-z_\.][A-z0-9\(\)_\.]+\s+)(.+)$/) {
        $line_modified = 1;
        for (@matrix_files) {
          print $_ $1;
        }
        $taxon_name = $1;
        my @tokens = $2 =~ /\s*\{[^\}]+\}|[^\{]/g;
        @tokens = map {$_ =~ s/\s+//; $_} @tokens; # Remove whitespace
        for (my $i = 0; $i < scalar(@tokens); $i++) {
          if ($tokens[$i] eq '-' && (substr $chartype[$i], 0, 1) =~ /([INT])/) {
            if ($1 eq 'N') { # Neomorphic character
              print $AMBIGUOUS  '?';
              print $AMBIGABS   '0';
              print $EXTRASTATE '0';
              print $INAPPLIC   '0';
            } elsif($1 eq 'I') { # Inversely-coded transformational character
              print $AMBIGUOUS  '?';
              print $AMBIGABS   '1';
              print $EXTRASTATE '1';
              print $INAPPLIC   '1';
            } elsif($1 eq 'T') { # Transformational character
              print $AMBIGUOUS  '?';
              print $AMBIGABS   '?';
              print $EXTRASTATE '9';
              print $INAPPLIC   '-';
            } else {
              warn "!! Unspecified character type at character $i\n";
              die;
            }
          } else {
            if ($tokens[$i] eq '9') {warn "!! Matrix already employs token '9'\n";}
            if ('{' eq substr $tokens[$i], 0, 1) {
              for (@matrix_files) {print $_ '?';} # Multi-state ambuguities not readily supported in R
            } else {
              for (@matrix_files) {print $_ $tokens[$i];}
            }
          }
        }
        for (@matrix_files) {print $_ "\n"};
        $taxon_name =~ s/^\s+//;
        $taxon_name =~ s/\s+$//;
        push @taxa, $taxon_name;
      }
      if (!$line_modified) {
        for (@matrix_files) {
          print $_ $line;
        }
      }
    }
    close EXTRASTATE;
    close AMBIGUOUS ;
    close AMBIGABS  ;
    close INAPPLIC  ;
    
    chdir $ROOT;
    # Now run analyses in TNT
    my @tnt_dirs = ("extraState/", "ambiguous/", "ambigAbsent/");
    foreach my $dir (@tnt_dirs) {
      my $tnt_file = $dir . $nexus_filename;
      system("tnt proc $tnt_file; run tnt_search.run $tnt_file;");
      
      open ($TNT_TREES, "<$tnt_file.tre") or warn ("!! Can't open $tnt_file.tre: $!\n");
      @tnt_trees = <$TNT_TREES>;
      close $TNT_TREES;
      pop @tnt_trees;
      shift @tnt_trees;
      
      open ($NEXUS_TREES, ">$tnt_file.nextrees") or warn ("!! Can't open $tnt_file.nextrees: $!\n");
      print $NEXUS_TREES "#NEXUS\nbegin taxa;\n\tdimensions ntax=" . scalar(@taxa) . ";\n\ttaxlabels";
      for (@taxa) {
        print $NEXUS_TREES "\n\t\t" . $_;
      }
      print $NEXUS_TREES "\n\t;\nend;\n";
      
      print $NEXUS_TREES "begin trees;\n";
      $i = 0;
      for (@tnt_trees) {
        ++$i;
        # Replace taxon numbers with taxon names
        s/([\s\(,])(\d+)(?= )/$1$taxa[$2],/g;
        s/, ?\)/)/g;
        # Annotations
        s/(=\S*)\//$1;/g;
        s/=(\S+)/[&Annot="$1"]/g;
        s/;/; /g;
        # Place commas between clades
        s/\)\(/),(/g;
        s/\]\s*\(/], (/g;
        # Separate multiple trees
        s/[\*; ]+;?$/;/;
        # Name trees
        s/^(.*?)([\d\w_]+)/tree tree$i = [&U] $1$2/;
        print $NEXUS_TREES "\t" . $_;
      }
      print $NEXUS_TREES "\nend;";
      close $NEXUS_TREES;
    }
  }
}
