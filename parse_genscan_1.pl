use lib "/home/hv270/perl5/lib/perl5";
use Bio::Tools::Genscan;
use Bio::SeqFeature::Gene::Transcript;
use Bio::Location::Simple;
use Bio::SeqFeatureI;

#arguments : slice_length, overlap_length,n_slices, outfile

if ($#ARGV != 3){
  print "Wrong Arguments\n";
  print "$#ARGV\n";
  foreach(@ARGV){
    print $_;
    print "\n";
  }
  exit;
}

my $outfile = @ARGV[3];
my $slice_length = @ARGV[0];
my $overlap_length = @ARGV[1];
my $n_genes = 0;

open (FILE, "> $outfile") || die "problem opening $outfile\n";

printf FILE 'gene_id,tag,gene_start,gene_end,strand,start,end,score';
print FILE "\n";

for ($i = 1;$i <= @ARGV[2];$i++){
  my $genscan = Bio::Tools::Genscan->new(-file => "/local/data/public/g1/Genscan/Dsim/slice_$i.txt");
  my $shift = ($i-1)*($slice_length - $overlap_length);
  while($gene = $genscan->next_prediction()) {
      # $gene is an instance of Bio::Tools::Prediction::Gene, which inherits
      # off Bio::SeqFeature::Gene::Transcript.
      @f= $gene->features();
      $loc= $gene->location();
      $n_genes++;
      foreach(@f){
        print FILE $n_genes;
        print FILE ",";
        print FILE $_->primary_tag();
        print FILE ",";
        print FILE $loc->start + $shift;
        print FILE ",";
        print FILE $loc->end + $shift;
        print FILE ",";
        print FILE $loc->strand;
        print FILE ",";
        print FILE $_->start;
        print FILE ",";
        print FILE $_->end;
        print FILE ",";
        print FILE $_->score();
        print FILE "\n";
      }
  }
  $genscan->close();
}

print "$n_genes genes\n";

close(FILE);
# essential if you gave a filename at initialization (otherwise the file
# will stay open)
