use File::Find;

$ROOT = "C:/Research/ExploreInapplicable/";
$MATRIX_DIR = "matrices/";
$CHAR_TYPE_DIR = "charType/";
print " Finding matrices in ". $ROOT . $MATRIX_DIR;
find (\&char_template, $ROOT . $MATRIX_DIR);
print "\nDone.";

sub char_template() {
	if (-f and /\.nex$/) {
    $infile = $_;
    $outfile_path = $ROOT . $CHAR_TYPE_DIR . $_;
    $outfile_path =~ s/\.nex$/.txt/i;
    $infile_path = $ROOT . $MATRIX_DIR . $infile;
    $i = 0;
    print "\n- Processing $infile ... ";
    if (-e $outfile_path) {
      print "$outfile_path already exists."
    } else {
      open (FILE, "<$infile_path") or warn " ERROR: Can't open $infile_path.\n";
      @lines = <FILE>;
      close FILE;
      
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
            $chars =~ s/\{[^\}]+?\}/?/g;
            $chars =~ s/\([^\)]+?\)/?/g;
            $chars =~ s/\[[^\]]+?\]/?/g;
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
      
      $i = 0;
      open (OUTFILE, ">$outfile_path") or warn "!! Can't open $outfile_path: $!\n";
      foreach $this_char (@char_is_inapp) {
        print OUTFILE ($this_char ? "NT" : "X")
            . ' [' . ++$i . "]\n";
      }
      close OUTFILE;
      print "Template written."
    }
  }
}
