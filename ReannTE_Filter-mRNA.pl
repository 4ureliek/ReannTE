#!/usr/bin/perl -w

#######################################################
# Author  :  Aurelie Kapusta
# email   :  4urelie.k@gmail.com
# version :  1.0
#######################################################
# Date    : v1.0, 30 Mar 2014
######################################################
use strict;
use warnings;
use Carp;
use Getopt::Long;
use Bio::SearchIO;
use Bio::SeqIO;

my $v = "1.0";
my $usage = "\nUsage, v$v:
     perl <scriptname.pl> -i <fa> [-b <blast-path>] [-e <XX>] [-forceB <X>] [-remote]
     OR
     perl <scriptname.pl> -i <fa> [-b <blast-path>] [-e <XX>] [-forceB <X>] [-db <fa>] [-dbt <XX>] [-bt <XX>]
     	
	(note that order of these options doesn't matter)
	(<...> denote places where user has to put its own file names / values)
	
	MANDATORY ARGUMENTS:
	 -i <fa>    => fasta file
	
	[OPTIONAL ARGUMENTS]:
	 -blast <path> => path = localisation of ncbi blast software
	                    if no path provided, path = /home/software/ncbi-blast-2.2.25+	 		            	          
	 -e <XX>       => XX = threshold, evalue (default = 10-10). It sets the minimum evalue to eliminate a sequence.
	 -forceB       => set x to chose how to behave if previous <fa>.blast.out exists
	                    x = 0 (default), chose this to avoid redoing the blast if <fa>.blast.out file already exists
	                    x = 1, chose this to save existing <fa>.blast.out (renamed), but still rerun blast
	                    x = 2, chose this to delete the pre-existing <fa>.blast.out file (therefore blast will be redone)
	 -remote       => use the -remote option of blast if you don't have the -db. This takes a while.
	 -db <fa>      => database to blast against [not relevant if -remote]
	 -dbt <XX>     => dbtype option of makeblastdb [default = nucl] [unrelevant if -remote]
	 -bt <XX>      => blast type [default = tblastx] [unrelevant if -remote]
	 	                  
	REQUIREMENTS:
	 - Blast software
	 - Bioperl

	WHAT IT DOES: 
	 Blast the (consensus) sequences 
	 against a db that can be defined only if not remote blast. If -remote, db = refseq_mrna
	 Then sequences that are filtered out  are only the ones without a class/family defined, or if this class or family are \"unclass\" or \"unknown\"
	 \n\n";
		 
#
# ---Variables and options
#
# declare variables for mandatory options
my $fa;
# declare variables + set default values for non mandatory options:
my $e = 0.0000000001; #10-10
my $blast_loc = "/home/software/ncbi-blast-2.2.25+";
my $forceB = 0;
my $remote = 0;
my $r;
my $db = "na";
my $dbtype = "nucl";
my $blast = "tblastx";
#get options
GetOptions ('i=s' => \$fa, 'db=s' => \$db, 'dbt=s' => \$dbtype, 'bt=s' => \$blast, 'b=s' => \$blast_loc, 'e=s' => \$e, 'forceB=s' => \$forceB, 'remote' => \$r);
$remote = 1 if $r;

# check that the fasta file is provided and is accessible
if ((! $fa) || (($remote == 0) && (! $db eq "na"))){
	die $usage;
}
if (! -e $fa){
	exit print "ERROR with input file, not provided or does not exists at location\n";
}
if (($remote == 1) && (! -e $fa)){
	exit print "ERROR check usage; remote not chosen + database defined does not exist at location\n";
}

#
# --- BLAST
#
my $out = "$fa.blast.out";
blast_on_split_seqs($fa,$blast_loc,$e,$db,$dbtype,$blast,$remote) unless ($forceB == 0 && -e $out);


#
# ---Filter hits
#
print " --- Filtering hits\n";
filter_mRNA($fa,$e);

print " --- Done\n\n";
exit;

##############################################################################################################
# --- SUBROUTINES
##############################################################################################################
#----------------------------------------------------------------------------
# split seqs
#----------------------------------------------------------------------------
sub blast_on_split_seqs {
	my ($fa,$blast_loc,$e,$db,$dbtype,$blast,$remote) = @_;	
	my $fa_obj = Bio::SeqIO->new(-file => $fa, -format => "fasta") or confess "Failed to open Bio::SeqIO object $fa $!\n";	
	
	#formatdb unless files exist or remote
	my $dbdb = "$db.nhr";
	$dbdb = "$db.phr" if ($dbtype =~ /prot/);
	unless ((-e "$dbdb") || ($remote == 1)) {
		my $fcmd = "$blast_loc/bin/makeblastdb -in $db -out $db -dbtype $dbtype";
		print "     Make blastdb, command line =\n     $fcmd\n";
		system $fcmd;
	}
	
	#now blast	
	if ($remote == 1) {
		my $fullout = "$fa.blast.out";
		my $dir;
		($fa =~ /\//)?($dir = "$fa.split"):($dir = "./$fa.split");
		system "rm -Rf $dir";
		mkdir $dir or confess "Failed to make dir $dir $!\n";
		print " --- split sequences and blast:\n";
		while( my $seq = $fa_obj->next_seq() ) {
			my $fullid = $seq->display_id;
			my $id = get_id($fullid);
			print "$id in progress\n";
			my $single = "$dir/$id";	
			my $single_obj = Bio::SeqIO->new(-file => ">$single", -format => "fasta") or confess "Failed to create Bio::SEQIO object $single $!\n";	
			$single_obj->write_seq($seq);
		
			#now the blast on the sequence
			my $out = "$single.blast.out";
			my $outprev = $out."prev";
			if ($forceB == 1 && -e $out) {
				system "mv $out $outprev";
				print "     $out exists, rename -> $outprev to save it before rerunning blast\n";
			}
			if ($forceB == 0 && -e $out) {
				print "     $out exists, skip blast\n";
			} else {
				my $cmd = "$blast_loc/bin/tblastx -db refseq_mrna -query $single -out $out -evalue $e -remote";
				print "     with command line =\n     $cmd\n";
				system $cmd;	
			}
		}
		#cat all outputs 	
		system "cat $dir/*blast.out > $fullout";
	} else {
		#no split
		my $out = "$fa.blast.out";
		my $cmd = "$blast_loc/bin/$blast -db $db -query $fa -out $out -evalue $e";
		print " --- blast with command line =\n     $cmd\n";
		system $cmd;	
	}	
}

#----------------------------------------------------------------------------
# get ID from fullname
#----------------------------------------------------------------------------
sub get_id {
	my $id = shift;
	$id =~ s/^(.*)#.*$/$1/ if ($id =~ /#/);
	return $id;
}

#----------------------------------------------------------------------------
# filter mRNA hits
#----------------------------------------------------------------------------
sub filter_mRNA {
	my ($fa,$e) = @_;

	my $out = "$fa.blast.out";
	open(my $out_fh, "<", $out) or confess "ERROR: could not open $out $!\n";	
	
	#parse blast output
	my ($reject,$retain) = parse_blast($out);
	print_fasta($fa,$reject,$retain);
}

#----------------------------------------------------------------------------
# filter mRNA hits
#----------------------------------------------------------------------------
sub parse_blast {
	my $blast = shift; #the blast result file
	my $out = "$blast.parsed.tab";
	my @reject = ();
	my @retain = ();
	my $blast_o = new Bio::SearchIO(-format => 'blast', -file => $blast); 
	open(my $out_fh,">",$out) or die "can not create file $out $!";
	print $out_fh 
	"#Qname\tQdescr\tHname\tHdescr\t%id\tevalue\tlength(including gaps)\tQlength\tHlength\tQstrand\tQstart\tQend\tHstrand\tHstart\tHend\n\n";
	while( my $result = $blast_o->next_result ) {
		while( my $hit = $result->next_hit ) {
			while( my $hsp = $hit->next_hsp ) {
				my $Qname = $result->query_name;
				my $id = get_id($Qname);
				my $Hdesc = $hit->description;		
				#push in the reject list, only if unclass repeat (or no class/fam set)
				push(@reject,$id) if (($Qname =~ /[Uu]nclass/) || ($Qname =~ /[Uu]nknown/) || ($Qname !~ /#/));
				#push in the retained list, if hit against transposase => these won't be eliminated
				push(@retain,$id) if ($Hdesc =~ /[Tt]ranspos/);
				# print to keep track
				print $out_fh 
					$result->query_name,"\t",
					$result->query_description,"\t",
					$hit->name,"\t",
					$hit->description,"\t",
					$hsp->percent_identity,"\t",
					$hsp->evalue,"\t",
					$hsp->length('total'),"\t",
					$hsp->length('query'),"\t",
					$hsp->length('hit'),"\t",
					$hsp->strand('query'),"\t",
					$hsp->start('query'),"\t",
					$hsp->end('query'),"\t",
					$hsp->strand('hit'),"\t",
					$hsp->start('hit'),"\t",
					$hsp->end('hit'),"\n"
				}  
		}
	}
	return (\@reject,\@retain);
}

#----------------------------------------------------------------------------
# print fasta files
#----------------------------------------------------------------------------
sub print_fasta {
	my ($fa,$list,$retain) = @_; #Here list = rejected stuff, retain = what to keep
	my @list = @{$list};
	my @retain = @{$retain};
	
	my $name = $fa;
	$name =~ s/\.fa$|\.fasta$//;
	
	#fa outputs
	my $keep = "$name.keep.fa";
	my $reject = "$name.reject.fa";
	
	my $fa_fa = Bio::SeqIO->new(-file => $fa, -format => "fasta") or confess "Failed to open Bio::SeqIO object $fa $!\n";	
	my $fa_keep = Bio::SeqIO->new(-file => ">$keep", -format => "fasta") or confess "Failed to create Bio::SeqIO object $keep $!\n";
	my $fa_reject = Bio::SeqIO->new(-file => ">$reject", -format => "fasta") or confess "Failed to create Bio::SEQIO object $reject $!\n";
	
	while( my $seq = $fa_fa->next_seq() ) {
		my $fullid = $seq->display_id;
		my $id = get_id($fullid);
		if ($id ~~ @list) { #if should be considered to be removed
			($id ~~ @retain)?($fa_keep->write_seq($seq)):($fa_reject->write_seq($seq)); #if is in retain list though, keep it
		}	
	}
}
