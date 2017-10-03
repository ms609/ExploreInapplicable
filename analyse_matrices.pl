use File::Find;

$ROOT = "C:/Research/ExploreInapplicable/";
$MATRIX_DIR = "matrices/";
$CHAR_TYPE_DIR = "charType/";
print " Finding matrices in ". $ROOT . $MATRIX_DIR;
find (\&char_template, $ROOT . $MATRIX_DIR);
print "\nDone.";

sub char_template() {
	if (-f and /V.*2008\.nex$/) {
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
      if (length ($_) && !($_ =~ /^[NTX]/i)) {warn ("!! Invalid chartype!");};
    }
    
    open (MATRIX, "<$raw_matrix_path") or warn " ERROR: Can't open $raw_matrix_path.\n";
    @nexus = <MATRIX>;
    close MATRIX;
    
    open ($AMBIGUOUS , ">$ambiguous_path"   ) or warn "!! Can't open $ambiguous_path: $!\n";
    open ($AMBIGABS  , ">$ambig_absent_path") or warn "!! Can't open $ambig_absent_path: $!\n";
    open ($EXTRASTATE, ">$extra_state_path" ) or warn "!! Can't open $extra_state_path: $!\n";
    open ($INAPPLIC  , ">$inapplicable_path") or warn "!! Can't open $inapplicable_path: $!\n";
    my @all_files = ($AMBIGUOUS, $AMBIGABS, $EXTRASTATE, $INAPPLIC);
    
    my $in_matrix = 0;
    foreach $line (@nexus) {   
      $line_modified = 0;
      if ($line =~ /^\s*matrix\s*$/i) {
        $in_matrix = 1;
      } elsif ($line =~ /;/) {
        $in_matrix = 0;
      } elsif ($in_matrix && $line =~ /^(\s*[A-z_]+\s+)(.+)$/) {
        $line_modified = 1;
        for (@all_files) {
          print $_ $1;
        }
        my @tokens = $2 =~ /\{[^\}]+\}|[^\{]/g;
        for (my $i = 0; $i < $#tokens; $i++) {
          if ($tokens[$i] eq '-' && (substr $chartype[$i], 0, 1) =~ /([NT])/) {
            if ($1 eq 'N') { # Neomorphic character
              print $AMBIGUOUS  '?';
              print $AMBIGABS   '0';
              print $EXTRASTATE '0';
              print $INAPPLIC   '0';
            } else { # Transformational character
              print $AMBIGUOUS  '?';
              print $AMBIGABS   '?';
              print $EXTRASTATE '9';
              print $INAPPLIC   '-';
            }
          } else {
            print $tokens[$i];
            for (@all_files) {print $_ $tokens[$i];}
          }
        }
        print "\n";
        for (@all_files) {print $_ "\n"};
      }
      if (!$line_modified) {
        for (@all_files) {
          print $_ $line;
        }
      }
    }
    die;
    
    while (!(shift (@lines) =~ /^\s*matrix\s*$/i)) {
      if (!$lines[0]) {
        warn " ERROR: No matrix line found";
        return;
      }
    }
    
    @char_is_inapp = ();
    $reading_matrix = 1;
    while ($reading_matrix) {
      $line = shift(@lines);
      if ($line =~ /;/) {
        $reading_matrix = 0;
      } else {
        if ($line =~ /^\s*[A-z_]+\s+(.+)$/) {
          $chars = $1;
          $chars =~ s/\{[\}]+\}/?/;
          $i = 0;
          for my $char (split //, $chars) {
            if ($char =~ /\-/) {
              $char_is_inapp[$i] = 1;
            } elsif ($char_is_inapp[$i] != 1) {
              $char_is_inapp[$i] = 0;
            }
            $i++;
          }
        }
      }      
    }
    close EXTRASTATE;
    close AMBIGUOUS ;
    close AMBIGABS  ;
    close INAPPLIC  ;
  }
}
