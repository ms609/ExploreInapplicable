use File::Find;

$ROOT = "C:/Research/ExploreInapplicable/";
$MATRIX_DIR = "matrices/";
$CHAR_TYPE_DIR = "charType/";
print " Finding matrices in ". $ROOT . $MATRIX_DIR;
find (\&char_template, $ROOT . $MATRIX_DIR);

print "\n Copying matrix for R";
open (NEXSRC, "<" . $dir . "/../lobo.nex") or warn " ERROR: Can't find NEXUS file $dir/../lobo.nex.\n";
@lines = <NEXSRC>;
close NEXSRC;
open (NEXR, ">$dir/../plot/lobo.nex")  or warn "!! Can't open ../plot/lobo.nex: $!\n";
for (@lines) {
  s/^([A-z_\s\(\)]{44})\s+/\1/g;
  s/\{\-,0\}/-/g;
  s/\{\-,1\}/-/g;
  s/\{0,2\}/A/g;
  s/\{3,4\}/B/g;
  s/\{1,2\}/C/g;
  s/\{0,1\}/D/g;
  print NEXR $_;
}
close NEXR;
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
      
      $i = 0;
      open (OUTFILE, ">$outfile_path") or warn "!! Can't open $outfile_path: $!\n";
      foreach $this_char (@char_is_inapp) {
        print OUTFILE ($this_char ? "NT" : "X")
            . ' [' . ++$i . "]\n";
      }
      close OUTFILE;
    }
  }
}
