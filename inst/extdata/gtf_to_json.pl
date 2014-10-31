
################################################################################
#
# read a GTF format file and write into json format
# The format of data structure is like:
# $gtf->{$gene_id}->{transcript}->{$transcript_id}->{exon}->{$exon_id}->{$start}
# $gtf->{$gene_id}->{transcript}->{$transcript_id}->{exon}->{$exon_id}->{$end}
# $gtf->{$gene_id}->{transcript}->{$transcript_id}->{exon}->{$exon_id}->{$i}
# $gtf->{$gene_id}->{transcript}->{$transcript_id}->{CDS}->{$exon_id}->{$start}
# $gtf->{$gene_id}->{transcript}->{$transcript_id}->{CDS}->{$exon_id}->{$end}
# $gtf->{$gene_id}->{transcript}->{$transcript_id}->{CDS}->{$exon_id}->{$i}
# $gtf->{$gene_id}->{transcript}->{$transcript_id}->{start}
# $gtf->{$gene_id}->{transcript}->{$transcript_id}->{end}
# $gtf->{$gene_id}->{transcript}->{$transcript_id}->{i}
# $gtf->{$gene_id}->{transcript}->{$transcript_id}->{type}
# $gtf->{$gene_id}->{strand}
# $gtf->{$gene_id}->{type}
# $gtf->{$gene_id}->{name}
# $gtf->{$gene_id}->{start}
# $gtf->{$gene_id}->{end}
# $gtf->{$gene_id}->{chr}
#
# Note: orders of genes are permutated because they are stored as hash.
#
# Usage:
#   perl gtf_to_json.pl gtf_file 1> gtf.json
#
# by Zuguang Gu
#
################################################################################


use JSON;
use strict;

my $gtf = read_gtf($ARGV[0]);

#print STDERR "writing to json...\n";
print to_json($gtf, { pretty => 1 } );

# note original order of genes or transcripts or exons will be disturbed.
sub read_gtf {
	my $file = shift;
				
	open my $fh, $file or die $!;
	my $gtf;
	my $exon_id;

	while(my $line = <$fh>) {
		next if($line =~/^#/);
		
		my @tmp = split "\t", $line;
		my $chr = $tmp[0];
		my $type = $tmp[2];
		my $start = $tmp[3];
		my $end = $tmp[4];
		my $strand = $tmp[6];
		
		my $gene_id = $tmp[8] =~/gene_id "(.*?)"/ ? $1 : die "cannot find gene_id!\n";
		my $transcript_id = $tmp[8] =~/transcript_id "(.*?)"/ ? $1 : $gene_id;
		my $transcript_type = $tmp[8] =~/transcript_type "(.*?)"/ ? $1 : "Unknown";
		my $gene_name = $tmp[8] =~/gene_name "(.*?)"/ ? $1 : $gene_id;
		my $gene_type = $tmp[8] =~/gene_type "(.*?)"/ ? $1 : "Unknown";
		
		# the new entry of a 'gene' is not defined by 'gene' column but by the description column
		if(!defined($gtf->{$gene_id})) {
			
			# create a new entry for this gene
			$gtf->{$gene_id} = {strand => "",
			                    name  => "",
								type  => "",
								start => 0,
								end => 0,
								chr => "",
								transcript => {},
								};
			$exon_id->{$gene_id} = {};
		}
		
		# only add new value if there is gene column
		if($type eq "gene") {
			$gtf->{$gene_id}->{strand} = $strand;
			$gtf->{$gene_id}->{name} = $gene_name;
			$gtf->{$gene_id}->{type} = $gene_type;
			$gtf->{$gene_id}->{start} = $start + 0;
			$gtf->{$gene_id}->{end} = $end + 0;
			$gtf->{$gene_id}->{chr} = $chr;
		}

		# because transcript_id in 'gene' is in fact gene_id
		# create an entry for transcript
		if($type ne "gene" and (! defined($gtf->{$gene_id}->{transcript}->{$transcript_id}))) {
			$gtf->{$gene_id}->{transcript}->{$transcript_id} = {exon => {},
			                                                    CDS => {},
																start => 0,
																end => 0,
																type => $transcript_type,
															};
			$exon_id->{$gene_id}->{$transcript_id} = 0;										
		}
	
		if($type eq "transcript") {
			$gtf->{$gene_id}->{transcript}->{$transcript_id}->{start} = $start + 0;
			$gtf->{$gene_id}->{transcript}->{$transcript_id}->{end} = $end + 0;
		}

		if($type eq "exon") {
			my @tt = split "; ", $tmp[8];
			$exon_id->{$gene_id}->{$transcript_id} ++;
			my $i_exon = $exon_id->{$gene_id}->{$transcript_id};
			$gtf->{$gene_id}->{transcript}->{$transcript_id}->{exon}->{"exon_$i_exon"} = { start => $start + 0,
			                                                                       end => $end + 0,
			                                                                       i => $exon_id->{$gene_id}->{$transcript_id} + 0,
			                                                                   };
		} elsif($type eq "CDS") {
			my @tt = split "; ", $tmp[8];
			my $i_exon = $exon_id->{$gene_id}->{$transcript_id};
			$gtf->{$gene_id}->{transcript}->{$transcript_id}->{CDS}->{"exon_$i_exon"} = { start => $start + 0,
			                                                                       end => $end + 0,
			                                                                      i =>$exon_id->{$gene_id}->{$transcript_id} + 0,
			                                                                  };
		}
		
		if($. % 100000 == 0) {
			#print STDERR "[serialize gtf] $. lines finished.\n";
		}
	}
	
	# test whether `gene` and `transcript` has start and end value...
	foreach my $gene_id (keys %$gtf) {
		my $min_transcript = 100000000000000;
		my $max_transcript = -1;
		foreach my $transcript_id (keys %{$gtf->{$gene_id}->{transcript}}) {
			my $min_exon = 100000000000000;
			my $max_exon = -1;

			foreach my $exon (keys %{$gtf->{$gene_id}->{transcript}->{$transcript_id}->{exon}}) {
				$min_exon = $min_exon < $gtf->{$gene_id}->{transcript}->{$transcript_id}->{exon}->{$exon}->{start} ? $min_exon : $gtf->{$gene_id}->{transcript}->{$transcript_id}->{exon}->{$exon}->{start};
				$max_exon = $max_exon > $gtf->{$gene_id}->{transcript}->{$transcript_id}->{exon}->{$exon}->{end} ? $max_exon : $gtf->{$gene_id}->{transcript}->{$transcript_id}->{exon}->{$exon}->{end};
			}
			
			if($gtf->{$gene_id}->{transcript}->{$transcript_id}->{start} == 0) {
				$gtf->{$gene_id}->{transcript}->{$transcript_id}->{start} = $min_exon + 0;
			}
			if($gtf->{$gene_id}->{transcript}->{$transcript_id}->{end} == 0) {
				$gtf->{$gene_id}->{transcript}->{$transcript_id}->{end} = $max_exon + 0;
			}

			$min_transcript = $min_transcript < $gtf->{$gene_id}->{transcript}->{$transcript_id}->{start} ? $min_transcript : $gtf->{$gene_id}->{transcript}->{$transcript_id}->{start};
			$max_transcript = $max_transcript > $gtf->{$gene_id}->{transcript}->{$transcript_id}->{end} ? $max_transcript : $gtf->{$gene_id}->{transcript}->{$transcript_id}->{end};
		}
		
		if($gtf->{$gene_id}->{start} == 0) {
			$gtf->{$gene_id}->{start} = $min_transcript + 0;
		}
		if($gtf->{$gene_id}->{end} == 0) {
			$gtf->{$gene_id}->{end} = $max_transcript + 0;
		}
	}

	return $gtf;
}
