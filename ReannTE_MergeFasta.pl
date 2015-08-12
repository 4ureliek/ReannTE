#!/usr/bin/perl -w

#######################################################
# Author  :  Aurelie Kapusta
# email   :  4urelie.k@gmail.com
#######################################################
use strict;
use warnings;
use Getopt::Long;
use Carp;
use Bio::DB::Fasta;
use Bio::SeqIO;
use Array::Unique;

my $version = "1.1";
my $changelog = "
#	- v1.0 = Mar 2014
#	- v1.1 = 06 Mar 2015
#             - remove log file, print STDERR instead
#             - deal with the -int added by Repeat Masker in Rname when _int
";
my $usage = "\nUsage, v$version:
     perl <scriptname.pl> -a <seqs_1.fa> -b <seqs_2.fa> [-p <x>] [-s <x>] [-forceRM <x>] [-gc <XX>] [-RM <path>] [-project <name>] [-CheckLow <XX>]
     	
	(note that order of these options doesn't matter)
	(<...> denote places where user has to put its own file names / values)
	
	MANDATORY ARGUMENTS:
	 -a <seqs_1.fa> => first fasta file
	 -b <seqs_2.fa> => second fasta file
	
	[OPTIONAL ARGUMENTS]:
	 -p <x>           => priority setting to favor or not one of the files when choice of sequence to keep
	                      x = a or b, give priority to file a or b when choice is not clear
	 		              x = no (default), both sequences will be kept
	 -s <XX>          => \"span\" corresponds to the minimum percentage of the sequence that is masked by another one to consider eliminating it
	                      The value [default = 80] will be used as a threshold to make choices on sequences to keep.
	                      For ex, if >XX% of sequenceA is masked by <XX% of sequenceB, sequenceB is kept. However, if <XX% of sequenceA is masked by <XX% of sequenceB, both are kept.
	 -forceRM        => set this to chose how to behave if previous .out exist
	                      x = 0 (default), chose this to avoid remasking if .out files already exist for files set as -a and -b
	                      x = 1, chose this to let RM check for existing .out (RM will move them if they do)
	                      x = 2, chose this to delete the pre-existing .out files (therefore masking will be redone)
	 -gc <XX>        => GC content (%) of the genome of the species considered, for use of good matrix in repeat masker               
	 -RM <path>      => path = localisation of repeat masker software
	 		              if no path provided, path = /home/software/RepeatMasker_405		            
	 -project <name> => name = will be in the name of the output files, including the merged fasta
	 		              if nothing provided, default = \"MergeFasta\"
	 -CheckLow <XX>  => chose this option to remove low complexity sequences before doing anything to merge libraries.
	 		              XX = threshold, in % (80% is advised). Set the minimum low complexity masked % required to eliminate the sequence.	 
	 -v              => verbose mode, make the script talks to you
	 -v              => print version if only option
	 -chlog          => print change log (updates)
	 -h|help         => Print this help
    		
	REQUIREMENTS:
	 - Repeat Masker software
	 - that ALL sequences have a unique name (e.g. name before the #)
	   if several different consensus have the same names between the 2 libraries this will create errors
	   you can use sed (see below) to add a number in front of all sequences of one of the files to avoid that issue in the case of merging 2 repclass outputs for ex.
	   sed 's/>/>1_/' seqs_1.fa > seqs_1.ok.fa 
	 
	WHAT IT DOES: 
	 This script merges two consensus libraries
	 - mask a with b (and b with a just to have access to it in case if needed)
	 - parses the masking outputs to evaluate overlaps
	 - make choices and flag sequences to keep or not. Note that all info are printed in an output => easy to go back and check, maybe adapt some parameters\n\n";
		 
#
# ---Variables and options
#
# declare variables for mandatory options
my $fa_a;
my $fa_b;

# declare variables + set default values for non mandatory options:
my $p = "no";
my $span = 80;
my $rm_loc = "/home/software/RepeatMasker_405";
my $forceRM = 0;
my $gc = "na";
my $outname = "MergeFasta";
my ($CheckLow,$help,$chlog,$v);
GetOptions ('a=s' => \$fa_a, 'b=s' => \$fa_b, 'p=s' => \$p, 's=s' => \$span, 'gc=s' => \$gc, 'RM=s' => \$rm_loc, 'forceRM=s' => \$forceRM, 'project=s' => \$outname, 'CheckLow=s' => \$CheckLow, 'chlog' => \$chlog, 'h' => \$help, 'help' => \$help, 'v' => \$v);

# check that the 2 fasta files are provided and they are accessible
die "\n version $version\n\n" if ((! $fa_a) && (! $fa_b) && (! $help)  && (! $chlog) && ($v));
die $usage if (((! $fa_a) || (! $fa_b)) || ($help));
die $changelog if ($chlog);
exit print STDERR "ERROR with input file, not provided or does not exists at location\n" if ((! -e $fa_a) || (! -e $fa_b));

# get path of lib_a => outputs there
my $path = path($fa_a);


#
# ---log
#
if ($v) {
	print STDERR "\n --- Script ReannTE_MergeFasta started, v$version:\n";
	print STDERR "     Low Complexity sequences will be eliminated, when sequences are more than $CheckLow % masked by low complexity\n" if ($CheckLow);
	print STDERR "     Files to be merged:\n";
	print STDERR "        - $fa_a\n";
	print STDERR "        - $fa_b\n";
	print STDERR "        -> Priority in choice given to file $p\n" if ($p ne "no");
	print STDERR "        -> No file will be favored (-p = $p)\n" if ($p eq "no");
	print STDERR "        -> min % of overlap (masking) = $span %\n";	
	print STDERR "     Repeat Masker software path = $rm_loc\n";
	print STDERR "        -> forceRM option set to $forceRM\n";
	print STDERR "\n";
}

#
# ---filter out low complexity if relevant
#
$fa_a = filter_low($fa_a,$rm_loc,$CheckLow,$v) if ($CheckLow);
$fa_b = filter_low($fa_b,$rm_loc,$CheckLow,$v) if ($CheckLow);


#
# ---mask the 2 libraries with each other
#
rm_RMout($fa_a,$fa_b) if ($forceRM == 2); #forceRM = 2, delete previous outputs (eg force remasking + do not bother with previous). $rm_check will be = 0
my $rm_check = 0;
my $out_a = $fa_a.".out";
my $out_b = $fa_b.".out";
$rm_check = 1 if ((-e $out_a) && (-e $out_b) && ($forceRM == 0)); #check for existence of previous output only if forceRM = 0

if (($rm_check == 0) || ($forceRM == 1)) {# => forceRM = 1, then just let RM do its stuff; if $forceRM ==1 + $rm_check = 1 means no previous outputs
	print STDERR " --- Repeat masking in progress (script is paused until completed)\n" if ($v);
	print STDERR "      -> with the following command line:\n" if ($v);
	my $ifgc;
	($gc eq "na")?($ifgc = ""):($ifgc = -gc $gc);
	my $cmda = "nohup $rm_loc/RepeatMasker $fa_a -lib $fa_b -e ncbi -xsmall -nolow $ifgc > $path/$outname.mask-a.nohup.log &";
	my $cmdb = "nohup $rm_loc/RepeatMasker $fa_b -lib $fa_a -e ncbi -xsmall -nolow $ifgc > $path/$outname.mask-b.nohup.log &";
	print STDERR "     $cmda\n" if ($v);
	print STDERR "     $cmdb\n" if ($v);
	system "$cmda";
	system "$cmdb";
	# pause script until maskings are both done
	sleep(5) until ((-e $out_a) && (-e $out_b)); #it was working without a time but then stopped working; find more elegant solution...??
	print STDERR "     Maskings done\n" if ($v);
	print STDERR "     Maskings done\n" if ($v);
} else {
	print STDERR " --- Repeat masking skipped, .out files exist in this folder\n" if ($rm_check == 1);
}


#
# ---now parse the repeat masker outputs to get information and "order" sequences
#
print STDERR " --- Parsing maskings\n" if ($v);
my ($complexity,$TEinfo,$mask_a,$mask_b,$data);

# ---get info (masked len for different repeats, nb of frgs, etc)
print STDERR "     Read .out outputs\n" if ($v);
print STDERR "      -> get more info for masking of sequences (perfragments + nb of fragments)\n" if ($v);
($mask_a,$data,$complexity,$TEinfo) = get_RMout($out_a,$data,$complexity,$TEinfo);
($mask_b,$data,$complexity,$TEinfo) = get_RMout($out_b,$data,$complexity,$TEinfo);

# ---get total length masked for each sequence THAT WAS IN .OUT
#    + get class fam etc + consensus lengths, key = Rname. Note that if some are not unique, it will be a problem.
print STDERR "     Read .masked outputs\n" if ($v);
print STDERR "      -> get info about consensus sequences: class, fam, length\n" if ($v);
print STDERR "      -> get total masked DNA (nt and %)\n" if ($v);
($complexity,$TEinfo) = get_masked_info_from_fa("$fa_a.masked",$complexity,$TEinfo,$v);
($complexity,$TEinfo) = get_masked_info_from_fa("$fa_b.masked",$complexity,$TEinfo,$v);


# ---get also complexity of the masking (complexity hash only here; both directions)
print STDERR "      -> get complexity, e.g. check how many repeats, family types, class types, are masking a given repeat\n" if ($v);
$complexity = get_complexity($mask_a,$data,$complexity,$TEinfo,$v);
$complexity = get_complexity($mask_b,$data,$complexity,$TEinfo,$v);


# ---loop in complexity to avoid undef values when unmasked sequences (e.g. non reciprocal maskings)
my @Ctypes = qw (desc_nb lst_Rname nb_Rname lst_Rfam nb_Rfam lst_Rclass nb_Rclass);
foreach my $seq (keys %{$complexity}) {
	foreach my $Ctype (@Ctypes) {
		$complexity->{$seq}{"$Ctype"} = 0 unless ($complexity->{$seq}{"$Ctype"});
	}	
}


# ---loop on all the 2 hashes to merge in one, just so printing output is correct => that reciprocal maskings won't appear twice.
print STDERR "     Merge maskings\n" if ($v);
my $mask = merge_mask($mask_a,$mask_b,$TEinfo,$complexity,$v);


# ---Then make choices between masked and masking one
print STDERR "     Make choices\n" if ($v);
my ($reject,$flag) = make_choice($mask,$complexity,$p,$span,$v);
 
#
# ---print outputfile 1: tabulated with all details
#
print STDERR " --- Print results 1/2\n" if ($v);
my $results = "$path/$outname.out.tab";
print_out($results,$mask,$complexity,$reject,$flag,$v);

my @list = ();
tie @list, 'Array::Unique';
foreach my $masked (keys %{$reject}) {
	foreach my $masking (keys %{$reject->{$masked}}) {
		push(@list,$reject->{$masked}{$masking}); #append ids (note they are not the full ids)
	}
}

#
# ---print outputfile 2: fasta file after decision (only the ones to keep + rejected ones)
#
print STDERR " --- Print results 2/2\n" if ($v);
my $fa_out = "$path/$outname.keep.fa";
my $fa_no = "$path/$outname.reject.fa";
print_fasta($fa_a,$fa_b,$fa_out,$fa_no,\@list,$v);

print STDERR "\n --- Script done\n\n" if ($v);
exit;




##############################################################################################################
# --- SUBROUTINES
##############################################################################################################
# SUMMARY OF THE SUBS
# - path
# - rm_RMout
# - get_id
# - get_class_fam
# - filter_low
#    - rw_fasta
# - get_RMout
#    - get_info_from_RMout_line
#    - get_TEinfo
# - get_masked_info_from_fa
# - merge mask
# - get_complexity
#    - get_homogeneity
# - make_choice
#    - priority_choice
# - print_out
######################################################
#----------------------------------------------------------------------------
# from a filename keep only its path - independently of the current dir
# my $path = path($fa_a);
#----------------------------------------------------------------------------
sub path {
	my($file) = shift;
	($file =~ /\//)?($file =~ s/(.*)\/.*$/$1/):($file = ".");
	return $file;
}

#----------------------------------------------------------------------------
# remove all outputfiles from repeat masker
# rm_RMout($fa_a,$fa_b)
#----------------------------------------------------------------------------
sub rm_RMout {
	my @files = @_;
	my @suffix = qw/cat log out ref tbl masked ori.out masked.index masked.infos/;
	foreach my $file (@files) {	
		foreach my $suffix (@suffix) {
			unlink "$file.$suffix" if (-e "$file.$suffix");
		}
	}	
}

#----------------------------------------------------------------------------
# move all outputfiles from repeat masker
# mv_RMout($fa,$folder)
#----------------------------------------------------------------------------
sub mv_RMout {
	my ($file,$folder) = @_;
	system "rm -Rf $folder" if (-e $folder);
	mkdir "$folder" or confess "\nERROR: Could not mkdir $folder $!\n";
	my @suffix = qw/cat log out ref tbl masked ori.out masked.index masked.infos/;
	foreach my $suffix (@suffix) {
		system "mv $file.$suffix $folder/" if (-e "$file.$suffix");
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
# get class,fam,classfam from fullname 
# (or not, but then no possible correction using Rname; unspecified => not the same thing.)
#----------------------------------------------------------------------------
sub get_class_fam {
	my $Rcf = shift;
	my $Rname = get_id($Rcf);
	$Rcf =~ s/^.*#(.*)$/$1/ if ($Rcf =~ /#/);
	$Rcf = $Rname if ($Rcf eq "Unspecified");
	my ($Rc,$Rf);
	if ($Rcf =~ /\//) {
		($Rc,$Rf) = split(/\//, $Rcf);
	} else {
		$Rc = $Rcf;
		$Rc=~ s/^(.*)\..*$/$1/;
		$Rf = $Rcf;
		$Rf =~ s/^.*\.(.*)$/$1/;
	}	
	return ($Rc,$Rf,$Rcf);
}

#----------------------------------------------------------------------------
# filter low complexity
# $fa_a = filter_low($fa_a,$rm_loc,$CheckLow,$v) if ($CheckLow);
# $fa_b = filter_low($fa_b,$rm_loc,$CheckLow,$v) if ($CheckLow);
#----------------------------------------------------------------------------
sub filter_low {
	my ($fa,$RM_loc,$p,$v) = @_;

	my $out = "$fa.out";
	print STDERR " --- Repeat masking low complexity in progress (script is paused until completed)\n" if ($v);
	print STDERR "      -> with the following command line:\n" if ($v);
	my $cmd = "nohup $rm_loc/RepeatMasker $fa -noint -e ncbi -xsmall > $fa.masklow.nohup.log &";
	print STDERR "         $cmd\n" if ($v);
	system ($cmd);		
	# pause until masking is done
	sleep(5) until (-e $out);
	
	#output
	my $name = $fa;
	$name =~ s/\.fa$|\.fasta$//;
	my $retained = "$name.nolow.fa";
	my $rejected = "$name.low.fa";
		
	if (-e "$fa.masked") {
		my $db = Bio::DB::Fasta->new("$fa.masked") or confess "ERROR: could not create Bio::DB::Fasta object from $fa.masked $!\n";	
		#create list of the ID of the file
		my @dbIDs = $db->get_all_ids();
	
		#fa output objects
		open(my $keep_fh, ">", $retained) or confess "ERROR: could not open to write $retained $!\n";	
		open(my $rejected_fh, ">", $rejected) or confess "ERROR: could not open to write $rejected $!\n";	

		print STDERR " --- Now filtering low complexity (using $p as % threshold)\n" if ($v);
		print STDERR "      => print kept sequences not low complexity in $retained\n" if ($v);
		print STDERR "      => print rejected low complexity sequences in $rejected\n" if ($v);		
		foreach my $Rfull (@dbIDs) {
			my $obj = $db->get_Seq_by_id($Rfull);		
			my $len = $obj->length;	#cons length	
			my $seq = $obj->seq; #sequence
			my $saveseq = $seq; #seq to print later
			$seq =~ s/[A-Z]//g; #remove anything NOT masked, including non conventional letters (just knowing that lowercase = masked)
			my $mlen = length $seq; #get masked length	
			my $mper = $mlen/$len*100;
			#print STDERR "$Rfull\t$mper\trejected\n" if ($mper > $p);
			#print STDERR "$Rfull\t$mper\tretained\n" if ($mper < $p);
			print $keep_fh ">$Rfull\n$saveseq\n" if ($mper < $p);
			print $rejected_fh ">$Rfull\n$saveseq\n" if ($mper > $p);
		}	
	} else {
		print STDERR " --- No low complexity or simple repeats in $fa ($fa.masked doesn't exist)\n" if ($v);
		$retained = "$fa";
	}
	rw_fasta($rejected);
	rw_fasta($retained);
	mv_RMout($fa,"$fa.RMlow"); #cleanup to avoid issue afterwards (other maskings)
	return $retained;
}

#----------------------------------------------------------------------------
# rewrite in clean fasta format, will replace previous file
#----------------------------------------------------------------------------
sub rw_fasta {
	my $fa = shift;
	my $fa_obj = Bio::SeqIO->new(-file => $fa, -format => "fasta") or confess "ERROR: failed to open Bio::SeqIO object $fa $!\n";

	my $realfa = "$fa.fasta";	
	my $realfa_obj = Bio::SeqIO->new(-file => ">$realfa", -format => "fasta") or confess "ERROR: failed to create Bio::SEQIO object $realfa $!\n";

	while( my $seq = $fa_obj->next_seq() ) {
		$realfa_obj->write_seq($seq);		
	}
	system "rm -f $fa";
	system "mv $realfa $fa"; 
}

#----------------------------------------------------------------------------
# open, read, RM .out file and store the infos, using 2 keys: masked ($Gname) and masking ($Rname).
# ($mask_a,$data,$complexity,$TEinfo) = get_RMout($out_a,$data,$complexity,$TEinfo);
# ($mask_b,$data,$complexity,$TEinfo) = get_RMout($out_b,$data,$complexity,$TEinfo);
#----------------------------------------------------------------------------
sub get_RMout {
	my ($RMout,$data,$complexity,$TEinfo) = @_;	
	my %mask = ();
	my %prev;
	open (RMOUT, "< $RMout") or confess "ERROR: Could not open $RMout $!\n";	
	LINE: while (<RMOUT>) {
		chomp (my $line = $_);	
		next LINE if (($line =~ /position|[sS]core|Gname/) || ($line !~ /\w/)); #skip headers and white lines
		my ($score,$div,$del,$ins,$Gfullname,$Gstart,$Gend,$Gleft,$Gmasked,$strand,$Rname,$Rclassfam,$Rstart,$Rend,$Rleft,$Rmasking,$ID) 
			= get_info_from_RMout_line($line); # => Rstart and Rleft, already corrected for strand)
		
		#Get TEinfo for both Gname and Rname unless it already exists (first filling of this hash)
		my $Rfullname = $Rname."#".$Rclassfam;
		my $Gname = get_id($Gfullname);
		$TEinfo = get_TEinfo($TEinfo,$Gfullname) unless ($TEinfo->{$Gname});
		$TEinfo = get_TEinfo($TEinfo,$Rfullname) unless ($TEinfo->{$Rname});
	
		# NOW PARSE
		# 1) Genomic stuff
		#Correct for overlapping maskings if relevant (overlap same masking repeat)
		my $overlap = 0;
		if ($prev{$Gname}{$Rname}{'G'} && ($prev{$Gname}{$Rname}{'G'} > $Gstart)) {
		 	$Gmasked = $Gend - $prev{$Gname}{$Rname}{'G'} +1;
		 	$overlap = 1;
		}
		$prev{$Gname}{$Rname}{'G'} = $Gend; #set for next line
		($mask{$Gname}{$Rname}{'Glen'})?($mask{$Gname}{$Rname}{'Glen'}+=$Gmasked):($mask{$Gname}{$Rname}{'Glen'}=$Gmasked);
		($mask{$Gname}{$Rname}{'nb'})?($mask{$Gname}{$Rname}{'nb'}++):($mask{$Gname}{$Rname}{'nb'}=1);		
		#keep all infos in case of several fragments of the same Rname
		my $i = $mask{$Gname}{$Rname}{'nb'};
		$data->{$Gname}{$Rname}{$i}="$div\t$ins\t$del\t$Gfullname\t$Gstart\t$Gend\t$strand\t$Rname\t$Rclassfam\t$Rstart\t$Rend\t$Rleft\t$ID";
	
		# 2) Consensus stuff
		$Rmasking-=$Gmasked if ($overlap == 1); #there was overlap in masking; this need also to be correcte in consensus masking length. Easier is to use that Gmasked amount, even if not the best.
		($mask{$Gname}{$Rname}{'Rlen'})?($mask{$Gname}{$Rname}{'Rlen'}+=$Rmasking):($mask{$Gname}{$Rname}{'Rlen'}=$Rmasking); 		
			
		#get total number of fragments
		($complexity->{$Gname}{'desc_nb'})?($complexity->{$Gname}{'desc_nb'}++):($complexity->{$Gname}{'desc_nb'}=1);		
	}
	return(\%mask,$data,$complexity,$TEinfo);
}

#----------------------------------------------------------------------------
# get values from RM output line, while (<RMOUT>)(=> split in subroutine)
#----------------------------------------------------------------------------
sub get_TEinfo {
	my ($TEinfo,$fullname) = @_;
	my $Rname = get_id($fullname);

	#Getting class and fam as well as possible + store them (masked sequence)
	my ($Rclass,$Rfam) = get_class_fam($fullname);	
	$TEinfo->{$Rname}{class} = $Rclass;
	$TEinfo->{$Rname}{fam} = $Rfam;
	$TEinfo->{$Rname}{fullname} = $Rname."#".$Rclass."/".$Rfam;
	
	return ($TEinfo);
}

#----------------------------------------------------------------------------
# get values from RM output line, while (<RMOUT>)(=> split in subroutine)
#----------------------------------------------------------------------------
sub get_info_from_RMout_line {
	my ($line) = @_;
	$line =~ s/^\s+//; #remove spaces in beginning of lines
	my ($score,$div,$del,$ins,$Gname,$Gstart,$Gend,$Gleft,$strand,$Rname,$Rclassfam,$Repstart,$Rend,$Repleft,$ID)= split(/\s+/,$line);
	
	#Get Rstart and Rleft, check strand
	my $Rleft = $Repleft;
	my $Rstart = $Repstart;
	$strand =~ s/C/-/;
	if ($strand eq "-") {
		$Rleft = $Repstart;
		$Rstart = $Repleft;
	}	
	my $Gmasked = ($Gend - $Gstart + 1);
	my $Rmasking = ($Rend - $Rstart + 1);
	
	#-int added at the end of some Rname. Pain in the butt since it means names don't match anymore...
	$Rname = $1 if ($Rname =~ /^(.*)-int$/);
	
	return ($score,$div,$del,$ins,$Gname,$Gstart,$Gend,$Gleft,$Gmasked,$strand,$Rname,$Rclassfam,$Rstart,$Rend,$Rleft,$Rmasking,$ID);
}

#----------------------------------------------------------------------------
# open, read, RM .masked file and store the infos
# ($complexity,$TEinfo) = get_masked_info("$fa_a.masked",$complexity,$TEinfo,$v);
#----------------------------------------------------------------------------
sub get_masked_info_from_fa {
	my ($fa,$complexity,$TEinfo,$v) = @_;
		
	my $infos = "$fa.infos";
	unless (-e $infos) {
		# index the genome and connect to the fasta file if relevant
		my $reindex;
		my $indexfile = "$fa.index";
		if (-e $indexfile) {
			$reindex = 0;
			print STDERR "        Fasta file previously indexed ($indexfile exists) - Skipping indexing step...\n" if ($v);
		} else {
			$reindex = 1;
			print STDERR "        Fasta file not indexed ($indexfile does not exists) - Indexing...\n" if ($v);
		}
		my $db = Bio::DB::Fasta->new($fa, -reindex=>$reindex) or confess "ERROR: could not create Bio::DB::Fasta object from $fa $!\n";
	
		#create list of the ID of the genome file
		my @dbIDs = $db->get_all_ids();
	
		#deal with getting seq lengths
		#first print (or not) fullname \t length \t mlength \t mper
		print STDERR "        -> lengths are being fetched + written in $infos\n" if ($v);	
		open (my $infos_fh, ">", $infos) or confess "ERROR: could not create file $infos $!\n\n";
		print $infos_fh "#Fullname\tseq_len\tmasked_len\tmasked_per\n";
		foreach my $Rfull (@dbIDs) {
			my $Rname = get_id($Rfull);
			if ($TEinfo->{$Rname}) { #do only for the ones that are relevant, e.g. masking or masked in both .out = are already keys in TEinfo
				my $obj = $db->get_Seq_by_id($Rfull);		
				my $len = $obj->length;	#cons length	
				my $seq = $obj->seq; #sequence
				$seq =~ s/[A-Z]//g; #remove anything NOT masked, including non conventional letters (just knowing that lowercase = masked)
				my $mlen = length $seq; #get masked length	
				my $mper = $mlen/$len*100;
				print $infos_fh "$Rfull\t$len\t$mlen\t$mper\n";			
				$TEinfo->{$Rname}{'len'}=$len;
				$complexity->{$Rname}{'desc_len'}=$mlen;
				$complexity->{$Rname}{'desc_per'}=$mper;
			}
		}		
	} else { #if file already there no need to index etc, just open the already printed info
		print STDERR "        -> lengths are being fetched from $infos (delete the file if infos should be fetched again from $fa)\n" if ($v);	
		open (my $infos_fh, "<", $infos) or confess "ERROR: could not open $infos $!\n";
		INFO: while (<$infos_fh>) {
			chomp (my $line = $_);
			next INFO if $line =~ /^#/;
			my ($Rfull,$len,$mlen,$mper) = split('\t',$line);
			my $Rname = get_id($Rfull);
			$TEinfo->{$Rname}{'len'}=$len;
			$complexity->{$Rname}{'desc_len'}=$mlen;
			$complexity->{$Rname}{'desc_per'}=$mper;		
		}
	}
	return ($complexity,$TEinfo);
}

#----------------------------------------------------------------------------
# merge mask - now that I have many infos, merge both maskings
# my $mask = merge_mask($mask_a,$mask_b,$TEinfo,$complexity,$v);
#----------------------------------------------------------------------------
sub merge_mask {
	my ($mask_a,$mask_b,$TEinfo,$complexity,$v) = @_;
	
	print STDERR "       = merge a and b mask hashes into one\n" if ($v);
	my $mask = ();
	# first, loop on mask_a 
	foreach my $masked (keys %{$mask_a}) {	
		foreach my $masking (keys %{$mask_a->{$masked}}) {			
			#store masked amount for $masked and masking amount for $masking		
			$mask_a->{$masked}{$masking}{'Gper'}=$mask_a->{$masked}{$masking}{'Glen'}/$TEinfo->{$masked}{'len'}*100;
			$mask_a->{$masked}{$masking}{'Rper'}=$mask_a->{$masked}{$masking}{'Rlen'}/$TEinfo->{$masking}{'len'}*100;
			
			if ($mask_b->{$masking}{$masked}) { #means that there is a reciprocal masking => check
				#check what complexity is higher => order in $mask
				if ($complexity->{$masked}{'desc_nb'} >= $complexity->{$masking}{'desc_nb'}) { #if these are =, keep seqs from file a first column
					$mask->{$masked}{$masking} = $mask_a->{$masked}{$masking};
					$mask->{$masked}{$masking}{'file'} = "a";
					$mask->{$masked}{$masking}{'reciprocal'}="yes";
				} else {
					$mask->{$masking}{$masked} = $mask_b->{$masking}{$masked};
					$mask->{$masking}{$masked}{'file'} = "b";
					$mask->{$masking}{$masked}{'reciprocal'}="yes";
				}				
			} else { #no reciprocal masking for this masked seq from a and masking from b -> mask_a kept
				$mask->{$masked}{$masking} = $mask_a->{$masked}{$masking};
				$mask->{$masked}{$masking}{'file'} = "a";
				$mask->{$masked}{$masking}{'reciprocal'} = "no";
			}
		}
	}
	# then deal with mask_b, but just sequences that were not in reciprocal masking
	foreach my $masked (keys %{$mask_b}) {
		foreach my $masking (keys %{$mask_b->{$masked}}) {
			#store masked amount for $masked and masking amount for $masking
			print STDERR "$masked, $masking\n" unless ($TEinfo->{$masking}{'len'});
			$mask_b->{$masked}{$masking}{'Gper'}=$mask_b->{$masked}{$masking}{'Glen'}/$TEinfo->{$masked}{'len'}*100;
			$mask_b->{$masked}{$masking}{'Rper'}=$mask_b->{$masked}{$masking}{'Rlen'}/$TEinfo->{$masking}{'len'}*100;
			
			unless ($mask_a->{$masking}{$masked}) {
				$mask->{$masked}{$masking} = $mask_b->{$masked}{$masking};
				$mask->{$masked}{$masking}{'file'} = "b";
				$mask->{$masked}{$masking}{'reciprocal'}="no";
			}	
		}
	}
	return $mask;
}

#----------------------------------------------------------------------------
# get percentage of sequence that is masked by other sequence
# $complexity = get_masked_details($mask,$data,$complexity,$TEinfo,$v);
#----------------------------------------------------------------------------
sub get_complexity {
	my ($mask,$data,$complexity,$TEinfo,$v) = @_;
	
	print STDERR "       - in progress\n" if ($v);	
	# foreach masked sequence:
	foreach my $masked (keys %{$mask}) {
		# check what is masking it and what is coverage
		foreach my $masking (keys %{$mask->{$masked}}) {
			my $currfam = $TEinfo->{$masking}{'fam'};
			my $currclass = $TEinfo->{$masking}{'class'};			
			$complexity = get_homogeneity($masking,"Rname",$complexity,$masked);
			$complexity = get_homogeneity($currfam,"Rfam",$complexity,$masked);
			$complexity = get_homogeneity($currclass,"Rclass",$complexity,$masked);
		}
	}
	return $complexity;	
}

#----------------------------------------------------------------------------
# check how many families, class etc
# $mask = get_homogeneity($currfam,"Rfam",$mask,$masked);
#----------------------------------------------------------------------------
sub get_homogeneity {
	my ($curr,$type,$complexity,$masked) = @_;
	my $t = "lst_".$type;
	my $n = "nb_".$type;
	if (($complexity->{$masked}{$t}) && ($complexity->{$masked}{$t} !~ /$curr/)) {		
		$complexity->{$masked}{$n}++;
		$complexity->{$masked}{$t} = $complexity->{$masked}{$t}.",".$curr;
	} else {
		$complexity->{$masked}{$n} = 1;
		$complexity->{$masked}{$t} = $curr;
	}
	return $complexity;
}

#----------------------------------------------------------------------------
# make choice
# my ($keep,$flag) = make_choice($mask,$complexity,$p,$span,$v);
#----------------------------------------------------------------------------
sub make_choice {
	my ($mask,$complexity,$p,$span,$v) = @_;

	print STDERR "       = compare %, length, complexity etc to make choices\n" if ($v);
	
	my $reject = ();
	my $flag = ();
	MASKED: foreach my $masked (keys %{$mask}) {
		MASKING: foreach my $masking (keys %{$mask->{$masked}}) {
			#set keep = "na" to see when choice was not made, to debug
			$reject->{$masked}{$masking}="na";
						
			#the 1:1 are easy case.
			if ($complexity->{$masked}{'nb_Rname'}==1 && $complexity->{$masking}{'nb_Rname'}==1) {
				# => check if they mask each other with > XX% (set by $span) consensus length	
				print STDERR "$masked\t$complexity->{$masked}{'desc_per'}%\t$masking\t$complexity->{$masking}{'desc_per'}%\n" if ($masked eq "1_R=1422");
				if ($complexity->{$masked}{'desc_per'}>$span && $complexity->{$masking}{'desc_per'}>$span){
					print STDERR "masking more than $span % => need to eliminate one\n" if ($masked eq "1_R=1422");
					#=> yes, so one needs to be eliminated
					$reject->{$masked}{$masking}=$masked if ($TEinfo->{$masked}{'len'} < $TEinfo->{$masking}{'len'});
					print STDERR "eliminated one will be $masked because ".$TEinfo->{$masked}{'len'}." < ".$TEinfo->{$masking}{'len'}."\n" if (($TEinfo->{$masked}{'len'} < $TEinfo->{$masking}{'len'}) && ($masked eq "1_R=1422"));
					if ($TEinfo->{$masked}{'len'} == $TEinfo->{$masking}{'len'}) {
						$reject = priority_choice($p,$mask,$masked,$masking,$reject);
					}
				} else {
					#If reciprocal masking
					#Note that sometimes non reciprocal doesn't mean it wouldn't be; if there are more than 2 sequences that are similar it will create issues...
					if ($mask->{$masked}{$masking}{'reciprocal'} eq "yes") {
						$reject->{$masked}{$masking}=$masking if ($complexity->{$masked}{'desc_per'} < $span && $complexity->{$masking}{'desc_per'} > $span); #masking is "included" in masked => discard masking one
						$reject->{$masked}{$masking}=$masked if ($complexity->{$masked}{'desc_per'} > $span && $complexity->{$masking}{'desc_per'} < $span); #masked is "included" in masking => discard masked one

						print STDERR "Reciprocal masking; eliminated one will be $masked because ".$complexity->{$masked}{'desc_per'}." > $span and ".$complexity->{$masking}{'desc_per'}." < $span\n" if (($complexity->{$masked}{'desc_per'} > $span) && ($complexity->{$masking}{'desc_per'}) < ($span) && ($masked eq "1_R=1422"));

						print STDERR "Reciprocal masking; eliminated one will be $masking because ".$complexity->{$masked}{'desc_per'}." < $span and ".$complexity->{$masking}{'desc_per'}." > $span\n" if (($complexity->{$masked}{'desc_per'} < $span) && ($complexity->{$masking}{'desc_per'} > $span) && ($masked eq "1_R=1422"));
					} else {
						# Deal with that	
											
					}	
				}
				next MASKED;
			}	
			#not 1:1.
			else {
				#check complexity of the masking sequences to help choice
				my $same_class = 0;
				my $same_fam = 0;
				$same_fam = 1 if (($complexity->{$masked}{'nb_Rfam'}==1) && ($complexity->{$masked}{'nb_Rname'}!=1));
				$same_class = 1 if (($complexity->{$masked}{'nb_Rclass'}==1) && ($complexity->{$masked}{'nb_Rname'}!=1));
				if ($complexity->{$masked}{'desc_per'}>$span) { # => case when masked one is >XX% and therefore will be excluded, unless priority is given to the file
					#Below, masking one is included (since there is another masking one contributing to the desc_per)
					$reject = priority_choice($p,$mask,$masked,$masking,$reject);
					#But this could suggest that masking sequences, even if kept, could be fragmented and be actually part of same TE. If true, complexity will be = 1 => flag
					#This is to consider only if repeats are included => >span				
					if (($reject->{$masked}{$masking} ne $masked) && ($complexity->{$masking}{'desc_per'}>$span)) { 
						$flag->{$masked}{$masking} = make_choice_flag($same_class,$same_fam);
					}
				} else { #not >XX% => masked is kept. Question = are the masking ones kept too? Could they be same repeat?
					if ($complexity->{$masking}{'desc_per'}>$span) { 
						$reject = priority_choice($p,$mask,$masked,$masking,$reject);
						
						#If retained, check that masking sequences could be fragmented and be actually part of same TE. If true, complexity will be = 1 => flag
						if ($reject->{$masked}{$masking} ne $masked) {
							$flag->{$masked}{$masking} = make_choice_flag($same_class,$same_fam);
						}
					}
				}
			}
		}
	}
	return ($reject,$flag);			
}

#----------------------------------------------------------------------------
# choice when equal "weight" for sequences, check if a file has priority or not
# 	using the fact that {$masked}{$masking}{'file'} countains info of what file this is from (since masking and masked would be inverted for file b)
# $keep = priority_choice($p,$mask,$masked,$masking,$keep);
#----------------------------------------------------------------------------
sub priority_choice {
	my ($p,$mask,$masked,$masking,$reject) = @_;
	if ($p eq "a") {
		($mask->{$masked}{$masking}{'file'} eq "b")?($reject->{$masked}{$masking}=$masked):($reject->{$masked}{$masking}=$masking);
	} elsif ($p eq "b") {
		($mask->{$masked}{$masking}{'file'} eq "a")?($reject->{$masked}{$masking}=$masked):($reject->{$masked}{$masking}=$masking);
	} else {
		$reject->{$masked}{$masking}="na";
	}
	return $reject;
}	

#----------------------------------------------------------------------------
# set flag if maybe repeats could be of the same element
# $keep = priority_choice($p,$mask,$masked,$masking,$keep);
#----------------------------------------------------------------------------
sub make_choice_flag {
	my ($same_class,$same_fam) = @_;
	my $flag;
	if (($same_fam == "1") && ($same_class == 1)){
		$flag = "Masking repeats = same family and class => unless repeats were manually curated consider that they belong to the same repeat";
	} elsif ($same_class == 1) {
		$flag = "Masking repeats = same class => even if not same family, unless repeats were manually consider that they belong to the same repeat";
	}
	return $flag;
}	

#----------------------------------------------------------------------------
# print output
# print_out($results,$mask,$complexity,$keep,$flag,$v);
#----------------------------------------------------------------------------
sub print_out {
	my ($out,$mask,$complexity,$reject,$flag,$v) = @_;
	
	print STDERR "       = print details in $out\n" if ($v);
	
	#open out file, print headers
	open(my $out_fh, ">", $out) or die "ERROR: could not create $out $!\n";	
	print $out_fh "#Note that masked sequences (first columns) will be from file defined as -a by default, unless several of these mask a single sequence from input file -b.\n\n";
	print $out_fh "#This_sequence\t.\t.\tis_masked_by:\t.\t.\t.This_seq_masks:\t.\tTotal_masked:\t.\tIf_reciprocal\tNumbers of masked:masking at various levels:\t.\t.\tSequence\n";
	print $out_fh "#Rname\tRclass\tRfam\tRname\tRclass\tRfam\t(nt)\t(%)\t(nt)\t(%)\tmasking:\tRclass\tRfam\tRname\tremoved\tComments\n\n";
	
	foreach my $masked (keys %{$mask}) {
		foreach my $masking (keys %{$mask->{$masked}}) {	
			#print name stuff
			print $out_fh "$masked\t$TEinfo->{$masked}{class}\t$TEinfo->{$masked}{fam}\t$masking\t$TEinfo->{$masking}{class}\t$TEinfo->{$masking}{fam}\t$mask->{$masked}{$masking}{'Glen'}\t$mask->{$masked}{$masking}{'Gper'}\t$complexity->{$masked}{desc_len}\t$complexity->{$masked}{desc_per}\t";
			#print if reciprocal
			($mask->{$masked}{$masking}{'reciprocal'})?(print $out_fh $mask->{$masked}{$masking}{'reciprocal'}."\t"):(print $out_fh "nd\t");
			
			#print homogeneity
			foreach my $key (sort keys %{$complexity->{$masked}}) {
				print $out_fh "$complexity->{$masked}{$key};$complexity->{$masking}{$key}\t" if $key =~ /^nb_/;
			}
			
			#print name of sequence that is kept
			($reject->{$masked}{$masking})?(print $out_fh $reject->{$masked}{$masking}."\t"):(print $out_fh "none\t");
			
			#print comments
			($flag->{$masked}{$masking})?(print $out_fh "$flag->{$masked}{$masking}\t"):(print $out_fh ".\t");
			
			#print line return
			print $out_fh "\n";
		}	
	}
}

#----------------------------------------------------------------------------
# print fasta files
# print_fasta($fa_a,$fa_b,$fa_out,$fa_no,\@list,$v);
#----------------------------------------------------------------------------
sub print_fasta {
	my ($fa_a,$fa_b,$fa,$fa_no,$list,$v) = @_;
	my @list = @{$list};
	
	print STDERR "       = print selected consensus sequences in $fa\n" if ($v);
	print STDERR "       = print rejected consensus sequences in $fa_no\n" if ($v);
	
	#cat fa_a and fa_b => open just one file
	my $fa_cat = "$fa_a.catb";
	system "cat $fa_a $fa_b > $fa_cat";
	
	my $fa_cat_obj = Bio::SeqIO->new(-file => $fa_cat, -format => "fasta") or confess "Failed to open Bio::SeqIO object $fa_cat $!\n";	
	my $fa_obj = Bio::SeqIO->new(-file => ">$fa", -format => "fasta") or confess "Failed to create Bio::SeqIO object $fa $!\n";
	my $fa_no_obj = Bio::SeqIO->new(-file => ">$fa_no", -format => "fasta") or confess "Failed to create Bio::SEQIO object $fa_no $!\n";
	
	while( my $seq = $fa_cat_obj->next_seq() ) {
		my $fullid = $seq->display_id;
		my $id = get_id($fullid);
		($id ~~ @list)?($fa_no_obj->write_seq($seq)):($fa_obj->write_seq($seq));		
	}
	unlink $fa_cat;
}


	