#!/usr/bin/perl -w

#######################################################
# Author  :  Aurelie Kapusta
# email   :  4urelie.k@gmail.com
# version :  1.0
#######################################################
# Date    : v1.0, 27 Mar 2014
######################################################
use strict;
use warnings;
use Carp;
use Getopt::Long;
use Bio::DB::Fasta;
use Bio::SeqIO;

my $v = "1.1";
my $usage = "\nUsage, v$v:
     perl <scriptname.pl> -i <fa> [-r <RMpath>] [-p <XX>]
     	
	(note that order of these options doesn't matter)
	(<...> denote places where user has to put its own file names / values)
	
	MANDATORY ARGUMENTS:
	 -i <fa>    => fasta file
	
	[OPTIONAL ARGUMENTS]:
	 -r <path>  => path = localisation of repeat masker software
	                  if no path provided, path = /home/software/RepeatMasker		          
	 -p <XX>    => XX = threshold, in % (default = 80%). It sets the minimum low complexity masked % required to eliminate the sequence
	                  
	REQUIREMENTS:
	 - Repeat Masker software, crossmatch engine
	 - Bioperl (Bio::DB::Fasta, Bio::SeqIO)

	WHAT IT DOES: 
	 This script masks low complexity (consensus) sequences in a fasta file
	 Eliminates the ones that are more than XX% masked (-p option)
	 2 fasta outputs: retained sequences and rejected sequences\n\n";
		 
#
# ---Variables and options
#
# declare variables for mandatory options
my $fa;

# declare variables + set default values for non mandatory options:
my $p = 80;
my $rm_loc = "/home/software/RepeatMasker";
GetOptions ('i=s' => \$fa, 'r=s' => \$rm_loc, 'p=s' => \$p);

# check that the fasta file is provided and is accessible
if (! $fa){
	die $usage;
}
if (! -e $fa){
	exit print "ERROR with input file, not provided or does not exists at location\n";
}

filter_low($fa,$rm_loc,$p);

print "\n --- DONE \n\n";
exit;

##############################################################################################################
# --- SUBROUTINES
##############################################################################################################
#----------------------------------------------------------------------------
# filter low complexity
#----------------------------------------------------------------------------
sub filter_low {
	my ($fa,$RM_loc,$p,$outname) = @_;

	my $out = "$fa.out";
	print "\n --- Repeat masking low complexity in progress, with command line:\n";
	my $cmd = "$rm_loc/RepeatMasker $fa -e crossmatch -noint -xsmall";
	print "     $cmd\n\n";
	
	print " --- start RM verbose -------------\n";
	system ($cmd);		
	print " --- end RM verbose ---------------\n";
	
	#output
	my $name = $fa;
	$name =~ s/\.fa$|\.fasta$//;
	my $keep = "$name.nolow.fa";
	my $rejected = "$name.low.fa";
		
	if (-e "$fa.masked") {
		my $db = Bio::DB::Fasta->new("$fa.masked") or confess "ERROR: could not create Bio::DB::Fasta object from $fa.masked $!\n";	
		#create list of the ID of the file
		my @dbIDs = $db->get_all_ids();
	
		#fa outputs		
		open(my $keep_fh, ">", $keep) or confess "ERROR: could not create $keep $!\n";			
		open(my $rejected_fh, ">", $rejected) or confess "ERROR: could not create $rejected $!\n";	
	
		#filter
		print "\n";
		print " --- Now filtering low complexity (using $p % threshold)\n";
		print "      => print kept sequences not low complexity in $keep\n";
		print "      => print rejected low complexity sequences in $rejected\n\n";
		print "Rfullname\tpercentage_low\tretained/rejected\n";
		foreach my $Rfull (@dbIDs) {
			my $obj = $db->get_Seq_by_id($Rfull);		
			my $len = $obj->length;	#cons length	
			my $seq = $obj->seq; #sequence
			my $saveseq = $seq; #seq to print later
			$seq =~ s/[A-Z]//g; #remove anything NOT masked, including non conventional letters (just knowing that lowercase = masked)
			my $mlen = length $seq; #get masked length	
			my $mper = $mlen/$len*100;
			print "$Rfull\t$mper\trejected\n" if ($mper > $p);
			print "$Rfull\t$mper\tretained\n" if ($mper < $p);
			print $keep_fh ">$Rfull\n$saveseq\n" if ($mper < $p);
			print $rejected_fh ">$Rfull\n$saveseq\n" if ($mper > $p);
		}	
	} else {
		print " --- No low complexity or simple repeats in $fa ($fa.masked doesn't exist)\n";
	}
	rw_fasta($keep);
	rw_fasta($rejected);
	mv_RMout($fa,"$name.RMlow"); #cleanup to avoid issue afterwards (other maskings)
}

#----------------------------------------------------------------------------
# rewrite in clean fasta format, will replace previous file
#----------------------------------------------------------------------------
sub rw_fasta {
	my $fa = shift;
	my $fa_obj = Bio::SeqIO->new(-file => $fa, -format => "fasta") or confess "Failed to open Bio::SeqIO object $fa $!\n";

	my $realfa = "$fa.fasta";	
	my $realfa_obj = Bio::SeqIO->new(-file => ">$realfa", -format => "fasta") or confess "Failed to create Bio::SEQIO object $realfa $!\n";

	while( my $seq = $fa_obj->next_seq() ) {
		$realfa_obj->write_seq($seq);		
	}
	system "rm -f $fa";
	system "mv $fa $fa"; 
}

#----------------------------------------------------------------------------
# move all outputfiles from repeat masker
# mv_RMout($fa,$folder)
#----------------------------------------------------------------------------
sub mv_RMout {
	my ($file,$folder) = @_;
	system "rm -Rf $folder" if (-e $folder);
	mkdir "$folder" or confess "\nERROR: Could not mkdir $folder $!\n";
	my @suffix = qw/cat log out ref tbl masked ori.out masked.index/;
	foreach my $suffix (@suffix) {
		system "mv $file.$suffix $folder/" if (-e "$file.$suffix");
	}	
}