use strict;
use warnings;

my$dir = '../../../VZV/ImmunoSEQ/VZVTCR_new/';

opendir(DH,$dir);
my@files = readdir(DH);
closedir(DH);

my%tcrseq;

$/ = "\r";

foreach my$file (@files){
	if($file =~ m/.+_\d+.txt/){
		open(IN,$dir.$file);
		print $file."\n";
		my%header;
		while(<IN>){
			my@l = split("\t");
			if($. == 1){
				for my$i (0 .. $#l){
					$header{$l[$i]} = $i;
				}
				print "Found amino acid at ".$header{"aminoAcid"}."\n";
				print "Found V at ".$header{"vGeneName"}."\n";
				print "Found J at ".$header{"jGeneName"}."\n";
			} else {
				if(not($l[$header{"aminoAcid"}] eq '')){
					if($l[$header{"aminoAcid"}] =~ m/\*/){next;}
					$tcrseq{join("\t",($l[$header{"vGeneName"}],$l[$header{"aminoAcid"}],$l[$header{"jGeneName"}]))} = 1;
				}
			}
		}
		
		close IN;
	}
}

open(OUT,">naive_VZVstudy.txt");
print OUT join("\n",keys(%tcrseq));
