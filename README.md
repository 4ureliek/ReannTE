ReannTE
=======
Last Update: 2014 08 06


Scripts to help transposable elements consensus sequences curation



========================================================
ReannTE_FilterLow.pl

     WHAT IT DOES: 
	 This script uses Repeat Masker to mask low complexity / simple repeats of the input fasta file
	 (for example, RepeatScout output)
	 
	 It eliminates the ones that are more than XX% masked (-p option)
	 2 fasta outputs: retained sequences and rejected sequences
	 	 
	 perl <scriptname.pl> -i <fa> [-r <RMpath>] [-p <XX>
	
     MANDATORY ARGUMENTS:
	 -i <fa>    => fasta file
	
     [OPTIONAL ARGUMENTS]:
	 -r <path>  => path = localisation of repeat masker software
	                  if no path provided, path = /home/software/RepeatMasker		          
	 -p <XX>    => XX = threshold, in % (default = 80%). It sets the minimum low complexity masked % required to eliminate the sequence
	                  
     REQUIREMENTS:
	 - Repeat Masker software, crossmatch engine
	 - Bioperl (Bio::DB::Fasta, Bio::SeqIO)

	