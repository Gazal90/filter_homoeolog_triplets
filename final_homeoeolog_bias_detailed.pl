#!usr/bin/perl
#
#Script to output details for homoeolog induction bias for each accepted triplet gene set under Fe stress
#USAGE: perl final_homeoeolog_bias_detailed.pl cuffdiffOK_file  acceptedHomoeologTRIPLETS_file OUTrejected_triplets > OUTrootFEvsControl_induction_bias

open(FE,$ARGV[0]);#-Fe cuffdiff file only OK
%category;
%ALL_ids;
while(<FE>)
{	chomp;
	if($_=~/Traes/)
	{	@geneDEinfo=split("\t",$_);chomp $geneDEinfo[2];chomp $geneDEinfo[13];chomp $geneDEinfo[9];
		@current_ids=split(",",$geneDEinfo[2]);
		foreach $id(@current_ids)
		{	if($_=~/Traes/)
			{	chomp $id; $a=$id;
				#print "$a\n";
			if(exists($ALL_ids{$a}))
			{$repeat++;	
				if($significance{$a}eq$geneDEinfo[13])
				{	$same++;#print "$ALL_ids{$a}\t$geneDEinfo[13]\n";
							#if((($foldchange{$a}<-1)&&($geneDEinfo[9]>1))||(($foldchange{$a}>1)&&($geneDEinfo[9]<-1)))
							#{print "$a\t$foldchange{$a}\t$geneDEinfo[9]\t$geneDEinfo[13]\n";}
				}
				else
				{	push @_reject,$a;
				}
			}
			else
			{	$significance{$a}=$geneDEinfo[13];	
				$ALL_ids{$a}=$geneDEinfo[12];
				$foldchange{$a}=$geneDEinfo[9];
			}
			}
		}
	}
}#	print @_reject;
#print %ALL_ids;
close FE;

@keys = ("nnn","yyy","ynn","nyn","nny","yyn","nyy","yny");
@values=("3flat","3up/3down","1up/1down","1up/1down","1up/1down","2up/2down","2up/2down","2up/2down");
@category{@keys}=@values;

open(TR,$ARGV[1]);#accepted homoeolog triplets
open(REJ,">",$ARGV[2]); #print rejected
while(<TR>)
{	chomp $_;@detail="";@values="";
	@triplet=split(" ",$_);
	$COUNT=0;
        $triplettype='';@exp='';
	$pos=0;$neg=0;
	foreach $gene (@triplet)
	{	if(exists($ALL_ids{$gene}))
		{	if($ALL_ids{$gene} < 0.05)
			{	$DE='y';
			}	else{$DE='n';}#print "$DE\n";
			#$triplettype=$triplettype.$DE;	
			if($DE eq 'y')
			{	#print "$DE\t";##########################
				if(($foldchange{$gene}) < -1)
				{	$neg++;
					push @detail, 'NEG';push @values,$foldchange{$gene};
				}
				elsif($foldchange{$gene} > 1)
				{	$pos++;
					push @detail, 'POS';push @values,$foldchange{$gene};
				}
				else{push @detail, 'NO';push @values,$foldchange{$gene};$DE='n';}
			}
			elsif($DE eq 'n'){push @detail, 'NO';push @values,$foldchange{$gene};}
			
			$triplettype=$triplettype.$DE;
			$COUNT++;
			foreach $rejectID(@_reject)
			{	chomp $rejectID; chomp $gene;
				if($gene ne $rejectID)
				{	#$COUNT++;print "$COUNT\n";
				}else{print REJ "DEFAULTER\t$_\n";
					$COUNT--;
					}
			}
		}
	}#print "###################################";
	if($COUNT==3)
	{	#print "$pos\t$neg\n";
		#print "$category{$triplettype}\t$_\n";#print "$triplettype\n";#print "$_\n";
		if($triplettype =~/y/ )
		{	if(($pos>0)&&($neg>0))
			{	print "OPPOSITE\t@detail\t@values\t$_\n";
			}
			else{	@specificCateg=split("/",$category{$triplettype});
				if($pos>0)
				{print "$specificCateg[0]\t@detail\t@values\t$_\n";
				}
				else{print "$specificCateg[1]\t@detail\t@values\t$_\n";}
			}
		}
		else{print "$category{$triplettype}\t@detail\t@values\t$_\n";}
	}
	else
	{#print REJ "$_\n";
	}
}
close TR;
$size=keys %ALL_ids;#print "$repeat\n$same";