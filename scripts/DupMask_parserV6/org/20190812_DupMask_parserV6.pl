#!/bin/env perl


#Zhaoshi Jiang
#add newHeader as ancestor reference file
#1-28-08 add function of HexTRgb color conversion function and chang color to Cytoband
#1-23-08 new format of the header of duplicons (SD123|chr1:123-456|p11)
#add size filter and range force to chain etc.
#9-24-07 new format of duplcicon file from Robert (MMU and PTR)
#7-23-07 find some dup dose not has anc information in the 92506subunitFamilyAnc
#force it to "NA 0 0" in the output file

#6-24-07 find some small pieces are due to Unmasked repeats, fitter small frags $opt_s

#5-25-07 add the function to chain frags with give ranges

#3-12-07 modified the code to fit whole genome data
#10-18-06 add -E when use -a one can only output dup,anc (6 columns).
#9-26-05 Tuesday add function to prepare show file for parasight
#9-25-06 add $opt_F (start coord file) for multiple fasta files
#9-25-06 incorporate liftover Anc result -u (subunit coord information bd35HaixuwBinary.4)
#9-21-06 add revese complement coord method
#8-15-06 Tuesday

#program description: Rober developet a dupMasker program to
#masker duplication sequences on various genome based on human
#dup consensus sequences library
#This code is to parse the output of that program and display
#the result by parasight



use Getopt::Std;

use List::Util 'min';
use List::Util 'max';

use vars qw ($opt_i $opt_o $opt_B $opt_c $opt_C $opt_w $opt_f $opt_F $opt_e  $opt_s $opt_S $opt_h $opt_E $opt_A $opt_R);

getopts("i:o:c:w:e:f:F:L:u:s:S:C:BAEhR");


if (!defined($opt_i) || defined($opt_h)) {

	PrintHelp();

}


#set default value#

$opt_o ||="$opt_i\.extra";
$opt_w ||=16;
$opt_e ||=65;
$opt_s ||=0; #minimal size filter  of dupmasker realignemt to be included
$opt_S ||=400000; #the frag size for each fasta 400k(hg17) or 500k (PTR2,MMU2)


# default color index files

if(defined ($opt_B)) { #colored by cytoband

$opt_c ||="/net/eichler/vol1/home/zhaoshi/hg17inform/cytoBandColor";

}


else { #colored by subunits

$opt_c ||="/net/eichler/vol1/home/zhaoshi/hg17inform/r.all.repeat.10K.nr.Color";


}

$opt_L ||="/net/eichler/vol1/home/zhaoshi/hg17inform/bd35ChrSize";
$opt_C ||=0; #range force to chain

####################################
#build a hash for default color list#
####################################

print "minimal size $opt_s, range force to chain $opt_C\n";


my %color=(
"chr1"=>"Misty Rose",
"chr2"=>"Black",
"chr3"=>"Dark Slate Gray",
"chr4"=>"Midnight Blue",
"chr5"=>"Cornflower Blue",
"chr6"=>"Medium Blue",
"chr7"=>"Deep Sky Blue",
"chr8"=>"Steel Blue",
"chr9"=>"Cyan",
"chr10"=>"Light Cyan",
"chr11"=>"Dark Green",
"chr12"=>"Spring Green",
"chr13"=>"Lawn Green",
"chr14"=>"Yellow",
"chr15"=>"Gold",
"chr16"=>"Rosy Brown",
"chr17"=>"Indian Red",
"chr18"=>"Brown",
"chr19"=>"Tomato",
"chr20"=>"Peru",
"chr21"=>"Red",
"chr22"=>"Hot Pink",
"chrX"=>"Deep Pink",
"chrY"=>"Maroon",
"chr_random"=>"Dark Violet",
"No hit"=>"Grey",
"ND"=>"Grey",
"NA"=>"Grey",
"Binary"=>"Grey",
"small"=>"Grey",
"No Net"=>"Grey"
);



#############################
######load input filed#######


use vars qw(%lenth %dupSize @input %anc %start_orient );


#input start and orient information

if ($opt_F) {load_data ("frag",$opt_F);} #load fragment start end and orient information


elsif ($opt_f) {

	my ($chr,$s,$e,$orient)=split(":",$opt_f);
	$start_orient{$chr}=join("\t",$s,$e,$orient);

}


#input ancestor and color information


if($opt_c){

		%color=(); #tosse the default colors
		load_data ("color",$opt_c); #input custmerized color


	} #or load color Index file for subunts

#input size information of chr and subunits

if ($opt_L) {load_data("Lenth",$opt_L);} #load the chromosome size information

#input segdupmasker result *.duplicons file

load_data ("input",$opt_i);


#chain the fragments within $opt_C

my $chainedRA=Chain(\@input);


#output different format output files

if ($opt_R) { #output RGB format file

	outputR($opt_o,$chainedRA);

}


else { #output Hexadecimal format

	output ($opt_o,$chainedRA);

}

###end of the main program###





##############################
##########subroutines########

sub outputR { #output color in RGB format ready for genome browser

	my $outfile=shift;

	my $rA=shift;

	my @fool=@$rA;

	open (OUT, ">$outfile") || die "Can not write to file $outfile\n";


	#print the header#

	print OUT "chr\tchrStart\tchrEnd\tSD\tscore\torient\tcolor\n";



	#end of print header


	foreach my $f (sort {$$a[0] cmp $$b[0] || $$a[1]<=>$$b[1] ||$$a[2]<=>$$b[2]} @fool) {

			my $query=join ("\t",$$f[0],$$f[1],$$f[2]);

			my $color="BEBEBE"; #grey color



			my $orient="";

			if ($$f[3] eq "F") {$orient="+";}

			elsif($$f[3] eq "R") {$orient="-";}




			my $subject=$$f[4];

			#new formate (SD123|chr1:123-456|p11)

			my $tmp=$$f[4];

			my ($subunit,$Anc,$AncS,$AncE,$cytoBand)=split(':|-|\|',$subject);

			$color=get_color("$Anc\t$AncS\t$AncE");

			#print "$color\n";

			$color=HexToRgb ($color); #covert to RGB color;

			#print "$color\n";


			#simplify the subject id

			$subject=~s/:0-0\|NA//; #simplify the NA

			if (!defined($opt_A)) {$subject=~s/\|.*//;}

			print OUT  "$query\t$subject\t0\t$orient\t$color\n";




		}


		close OUT;



}



sub output { #convert small frags to uniform color (black)


	my $outfile=shift;

	my $rA=shift;

	my @fool=@$rA;

	open (OUT, ">$outfile") || die "Can not write to file $outfile\n";


	#print the header#

	if ($opt_E) {

		print OUT "chr\tchrStart\tchrEnd\torient\tRepeat\tcolor\twidth\toffset\n";
	}

	else {

		print OUT "qChr\tqStart\tqEnd\tOrient\tSubunit\taChr\taStart\taEnd\n";
	}

	#end of print header


	foreach my $f (sort {$$a[0] cmp $$b[0] || $$a[1]<=>$$b[1] ||$$a[2]<=>$$b[2]} @fool) {

			my $query=join ("\t",$$f[0],$$f[1],$$f[2]);

			my $color="BEBEBE"; #grey color

			my $orient=$$f[3];

			my $subject=$$f[4];


			my $size=$$f[2]-$$f[1]+1;

			if ($size<$opt_s && $opt_s>0) { #make small dup uniform look and rm ID

			     $color="#000000"; #black

			     $orient="";

			     $subject="";


			}

			else { #large dup



				#new formate (SD123|chr1:123-456|p11)

				my $tmp=$$f[4];

				my ($subunit,$Anc,$AncS,$AncE,$cytoBand)=split(':|-|\|',$subject);

				$color=get_color("$Anc\t$AncS\t$AncE");

			}


				#simplify the subject id

				$subject=~s/:0-0\|NA//; #simplify the NA

				if (!defined($opt_A)) {$subject=~s/\|.*//;}


			#print out differen result
			if($opt_E) {

				print OUT  "$query\t$orient\t$subject\t$color\t$opt_w\t$opt_e\n"


			} #$opt_E

			else {

				print OUT "$query\t$orient\t$subject\n";

			}


		}


		close OUT;



}


sub Chain { #force to chain of subunit fragements within give range default 10kb

	my $rA=shift;
	my @fool=@$rA;
	my @result=();
	my $qname="";
	my $s=0;
	my $e=0;
	my $orient="";
	my $subunit="";
	my $size=0;

	foreach my $f (sort{$$a[0] cmp $$b[0] || $$a[1] <=> $$b[1] || $$a[2]<=>$$b[2] || $$a[4] cmp $$b[4]} @fool) {


			#print "$$f[0]\n";

			if ($qname eq "") {

						$qname=$$f[0];
						$s=$$f[1];
						$e=$$f[2];
						$orient=$$f[3];
						$subunit=$$f[4];
			} #first record


			elsif ($qname eq $$f[0] && min($e,$$f[2])+$opt_C>=max($s,$$f[1]) && $subunit eq "$$f[4]" && $orient eq "$$f[3]") {
						$s=min($s,$$f[1]);
						$e=max($e,$$f[2]);

			} #elsif


			else { #new subunit or out of range


				push @result, [$qname,$s,$e,$orient,$subunit];
				$qname=$$f[0];
				$s=$$f[1];
				$e=$$f[2];
				$orient=$$f[3];
				$subunit=$$f[4];


		} #else


	} #my $f


		push @result, [$qname,$s,$e,$orient,$subunit];



	return \@result;


}



sub prepare_showfile {

	open (SHOW,">showfile")|| die "Can not write to file showfie\n";

	foreach my $key (sort keys %start_orient) {


		my ($s,$e,$orient)=split("\t",$start_orient{$key});

		#correct orient
		if ($s>$e) {my $sth=$e; $e=$s;$s=$sth;}

		my $key2=$key;

		$key2=~s/_[0-9]+$//;

		print SHOW "$key\t$lenth{$key2}\t$s\t$e\n";



	} #foreach



} #end of prepare

sub get_color { #3-13-07

 	my $fool=shift;

 	#print "$fool\n";

 	my $color="#BEBEBE"; #grey

 	my ($anc,$s,$e)=split("\t",$fool);

 	#if defined $opt_c  (chr,s,e,color)

 	if (defined($opt_c)) { #with ancestor file and color index file

		foreach my $key (sort keys %color) {

		 	#die "$key\t$color{$key}\n";

		 	my ($runChr,$runS,$runE)=split("\t",$key);


		 	if ($anc eq $runChr && min($e,$runE)>max($s,$runS)) {

		 		#die "$anc $s $runChr $runS";
		 		$color=$color{$key};
		 		last;

		 	}



		 } #foreach


 	} #if


 	else { #no color index files assign the color based on ancestral chr

 			if($anc=~/_random/) {$color=$color{"chr_random"};}

 			elsif(exists $color{$anc}){$color=$color{$anc};}

 	}



 	return ($color);




}

sub load_data {

	my $index=shift;
	my $file=shift;

	open (IN,"$file") || die "Can not read the file $file\n";

	while (<IN>) {

		s/\r\n/\n/;

		chomp;

		my @tmp=split('\s',$_);

		if ($index eq "Lenth") { #input chromosome Length file

			next if ($tmp[1]=~/\D/); #skip the header
			$lenth{$tmp[0]}=$tmp[1];

		} #end of if

		elsif ($index eq "lenth") { #input dup subunit size information

			next if ($tmp[1]=~/\D/); #skip the header
			$dupSize{$tmp[0]}=$tmp[1];

		}


		elsif ($index eq "frag") { #start coord and orient

			$start_orient{$tmp[0]}=join ("\t",$tmp[1],$tmp[2],$tmp[3]);

		}

		elsif ($index eq "color") {

			my $key="";


			if (scalar(@tmp)==2) { #tmp_S* #A168B9

				$key=$tmp[0];

				$color{$key}=$tmp[1];
			}

			elsif (scalar(@tmp)==4){ #chr1    465     167280  #A168B9

				$key=join("\t",$tmp[0],$tmp[1],$tmp[2]);

				#print "$key\t$tmp[3]\n";

				$color{$key}=$tmp[3];

			}


			elsif (scalar(@tmp)==5) { #chr1    0       2299999 chr1_p36.33     #B88444

				$key=join("\t",$tmp[0],$tmp[1],$tmp[2]);

				#print "$key\t$tmp[3]\n";

				$color{$key}=$tmp[4];


			}



		}

		elsif ($index eq "input") {

				my $chr="";
				my $subOrient="F";
				my $line=$_;
				$line=~s/\r\n/\n/;
				chomp ($line);

				#254 18.07 1.20 2.44 chr2dna 28197 28279 (561442) C tmp-1_s11198 (1) 2012 1931 7
				#6267 9.83 3.39 0.33 AC097264.4 26344 27228 (146492) tmp-823_s2 81 992 (25) 1

				#whole genome data (Hg17)
				#55618 0.42 0.00 0.13 chr22_036 30002 36002 (363998) tmp-1_s6670 2 5994 (0) 1

				#(PTR2)
				#12009 7.87 0.19 0.64 chr1_193 2 1563 (458521) tmp-1_s3635 660 2214 (1) 1


				#(MMU2)

				#1078 33.50 4.09 1.62 split368 207450 208231 (291769) C tmp-1_s17347 (30438) 97511 96711 13


				# MMU2 (v5)
				#1979 13.04 0.67 0.67 chr8.fa_314 79360 79658 (320342) C tmp-1_s18713 (6698) 6394 6096

				#new format of duplicon header: (SD123|chr1:123-456|p11)

				if ($line=~m/ C /) { #reverse complement

					$subOrient="R";
					$line=~s/ C / /g;

				}


				$line=~s/^[\s]+//;

				my @tmp=split(/[\s]+/,$line);


				#whole genome data
				#55618 0.42 0.00 0.13 chr22_036 30002 36002 (363998) tmp-1_s6670 2 5994 (0) 1

				if($tmp[4]=~/^chr[0-9XYrandomU\.fa]+_/) {

					$tmp[4]=~s/\.fa//; #rm the fa 1-22-08

					#print "$tmp[4]\n";

					($chr,$trueStart,$trueEnd,$trueOrient)=true_start_orient2($tmp[4],$tmp[5],$tmp[6],$subOrient); #3-12-07 add to fit whole genoem data

				}

				#1078 33.50 4.09 1.62 split368 207450 208231 (291769) C tmp-1_s17347 (30438) 97511 96711 13

				elsif ($tmp[4]=~/split/) { #new input format 9-20-07 for 500kb frags


					#chr1.duplicon
					#get chr name
					my $tmp="$file";
					my @sth=split('\.',$tmp);
					$chr="chr$sth[0]";
					#get start end orient and subunit

					($trueStart,$trueEnd,$trueOrient)=true_start_orient3($tmp[4],$tmp[5],$tmp[6],$subOrient);


				}


				else { #clone data

					($chr,$trueStart,$trueEnd,$trueOrient)=true_start_orient($tmp[4],$tmp[5],$tmp[6],$subOrient,\%start_orient);

				}


				# "$chr,$trueStart,$trueEnd,$trueOrient,$tmp[8]\n";

				push @input,[$chr,$trueStart,$trueEnd,$trueOrient,$tmp[8]];


		} #end of elsif



	} #while


} #end of load_data

sub true_start_orient3 { #9-20-07 for new format of input 500kb frags from MMU and PTR

	my $chr_fragN=shift;
	my $sub_start=shift;
	my $sub_end=shift;
	my $sub_orient=shift; #the subunit orient
	#1078 33.50 4.09 1.62 split368 207450 208231 (291769) C tmp-1_s17347 (30438) 97511 96711 13

	$chr_fragN=~s/split//;

	my $frag_s=$chr_fragN*$opt_S; #chr fasta was fragmentated into 400kb frag
	my $frag_e=$chr_fragN*$opt_S;


	my $trueStart=$frag_s+$sub_start;
	my $trueEnd=$frag_s+$sub_end;
	my $trueOrient=$sub_orient;

	#die "$chr,$trueStart,$trueEnd,$trueOrient\n";

	return ($trueStart,$trueEnd,$trueOrient);



}


sub true_start_orient2 { #3-12-07 for whole genome data chr*_001


	my $chr_fragN=shift;
	my $sub_start=shift;
	my $sub_end=shift;
	my $sub_orient=shift; #the subunit orient

	#chr15_045
	#chr15_random_000
	#chrY_50 (PTR2)
	my $chr="";
	my $fragN=0;

	if ($chr_fragN=~/random/) {

		$chr=substr ($chr_fragN, 0,-4);
		$fragN=substr ($chr_fragN, -3);


	}

	else {

		($chr,$fragN)=split("_",$chr_fragN);


	}

	#print "$chr $fragN $opt_S\n";



	my $frag_s=$fragN*$opt_S; #chr fasta was fragmentated into 400kb frag
	my $frag_e=$fragN*$opt_S;


	my $trueStart=$frag_s+$sub_start;
        my $trueEnd=$frag_s+$sub_end;
	my $trueOrient=$sub_orient;
	return ($chr,$trueStart,$trueEnd,$trueOrient);





}


sub true_start_orient {

	my $chr=shift;
	my $sub_start=shift;
	my $sub_end=shift;
	my $sub_orient=shift; #the subunit orient
	my $rH=shift;
	my %fool=%$rH;

	my ($frag_s,$frag_e,$frag_orient)=(0,0,"F"); #fragment start end and orient

	my ($trueStart,$trueEnd,$trueOrient)=(0,0,"F");

	if (exists $fool{$chr}) {

		($frag_s,$frag_e,$frag_orient)=split("\t",$fool{$chr});


	}


	if ($frag_orient eq "F") { #fragment in forward orientation

		$trueStart=$frag_s+$sub_start;
		$trueEnd=$frag_s+$sub_end;
		$trueOrient=$sub_orient;

	} elsif ($frag_orient eq "R"){

		$trueStart=$frag_e-$sub_start;
                $trueEnd=$frag_e-$sub_end;
		if($sub_orient eq "F") {

			$trueOrient="R";

		} elsif ($sub_orient eq "R"){

			$trueOrient="F";
		}

	}

	return ($chr,$trueStart,$trueEnd,$trueOrient);

}


####color coversion HextoRGB function start here####

sub HexToRgb {


	##fffafa;

	my $hex=shift;
	my $r=TranColor(substr($hex,1,1))*16+TranColor(substr($hex,2,1));
	my $g=TranColor(substr($hex,3,1))*16+TranColor(substr($hex,4,1));
	my $b=TranColor(substr($hex,5,1))*16+TranColor(substr($hex,6,1));

	my $rgb=join(",",$r,$g,$b);

	return ($rgb);

}




sub TranColor {

	my $color=shift;

 	if ($color eq "A") {$color=10;}

 	elsif ($color eq "B") {$color=11;}

 	elsif ($color eq "C") {$color=12;}

 	elsif ($color eq "D") {$color=13;}

 	elsif ($color eq "E") {$color=14;}

 	elsif ($color eq "F") {$color=15;}

 	return $color;


 }


#######color funtion end of here####



sub PrintHelp {

my $program=$0;
$program=~s/.*\///;
print "$program parse the segDupmask (Version6) output file (new format)
-i [input from dupmasker *.duplicons]
-c [color Index file, 2 colums:binary_s0 #B55531 or 4 columns:anc,s,e #B55431]
-f [only one frag, start,end,orient: chr:s:e:F/R]
-C [range force to chain,default no chaining]
-F [a list of frags start,end,orient information, default 0:0:F]
-L [chromosome Length file]
-E [flag to output extra file, when use -a But NOT -E,one only output dup,anc]
-w [width of bar, defaut 16]
-e [extra file start position,default 65]
-s [minimal size cutoff of dupmaker hits to be incldued, default 0bp]
-S [fasta fragment size defautl 400kb(HG17), 500kb for PTR2/MMU2]
-A [switch to output Ancestor information too, default only subunit ID]
-R [swithc to covert color in RGB format ouput file ready for browser]
-B [color the duplicon based cytoBand location,default by subunits]
-o [output file]
-h [help file]\n";
exit;

}


