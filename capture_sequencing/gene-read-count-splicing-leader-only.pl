#!/usr/bin/perl -w

#2020-05-07

#1) only count reads with covid-19 leader seq 
#(The Architecture of SARS-CoV-2 Transcriptome: splicing 5' start at: 55 and 85 bp in the paper
# we used slight loose leader range for splicing 5': 34-85 based on the leader sequence
#2) gene read count: 3' of splice read fall upstream
#                    3' of splice read inside the gene, then count as next gene
#3) also print splice juntion coordinates in SJ.txt: ref  5'-splice  3'-splice  count-reads

#only read CIGAR string and MD tags
#1 QNAME Query (pair) NAME 
#2 FLAG bitwise FLAG 
#3 RNAME Reference sequence NAME 
#4 POS 1-based left most POSition/coordinate of clipped sequence 
#5 MAPQ MAPping Quality (Phred-scaled) 
#6 CIAGR extended CIGAR string
#      M Alignment match (can be a sequence match or mismatch)
#      I Insertion to the reference
#      D Deletion from the reference
#      N Skipped region from the reference
#      S Soft clip on the read (clipped sequence present in <seq>)
#      H Hard clip on the read (clipped sequence NOT present in <seq>)
#      P Padding (silent deletion from the padded reference sequence)

#7 MRNM Mate Reference sequence NaMe ('=' if same as RNAME) 
#8 MPOS 1-based Mate POSistion 
#9 ISIZE Inferred insert SIZE 
#10 SEQ query SEQuence on the same strand as the reference 
#11 QUAL query QUALity (ASCII-33 gives the Phred base quality) 
#12 OPT variable OPTional fields in the format TAG:VTYPE:VALUE 

use strict;
my($sample, $in,$gtf, $out,$min_splice) =@ARGV;
if(!$sample||!$in||! $gtf||! $out ||$in eq $out||$gtf eq $out){
    die "\tUsage: sample-name  input-sam/bam(sam/bam)  gene-gtf  out-gene-read-count  min-splice-size(default 20)\n";
}

if(-f $out){
    #die "Out file exists (please remove it if you want to run this script): $out\n";
}

if(!$min_splice){
    $min_splice = 20;
}
my ($leader_start, $leader_end) =(34,85);  #5' splicing start in leader

my $out_sj = $out;
$out_sj =~ s/\.txt$//;
$out_sj =~ s/\.tab$//;
$out_sj = $out_sj.".SJ.txt";                #splice junction file: ref 5'-splice 3'-splice count-reads

#parse gtf to get gene plus upstream coordinates
my %HOGU = ();           #gene plus upstream coordinates
my @genes =();

open(GTF, $gtf) or die "can't open gene GTF file: $gtf\n";
my($cnt_gene) = (0);
my $pre_end = '';
while(<GTF>){
    next if (/^\#/);
    if(/\S/){
	chomp;
	my($id,$source,$type,$s,$e,$src,$strand,$phase,$attribute) = split '\t', $_;
	if($type eq 'gene'){
	    if($attribute !~ /(\"5\'UTR\")|(\"3\'UTR\")/){    #hardcode for covid-19
		my ($gene) = $attribute =~ /gene ([^;]+);/;  $cnt_gene++;
		if(!$gene){
		    die "unknown gft format (no gene name as gene NAME;): $_\n";
		}
		my $upstream = $s;
		if($pre_end){
		    $upstream  = $pre_end +1;
		}else{
		    $upstream = 1;   #first gene, upstream =  1
		}
		$HOGU{$gene}{'upstream'} =  $upstream;
		$HOGU{$gene}{'s'} =  $s;
		$HOGU{$gene}{'e'} =  $e;
		$HOGU{$gene}{'chr'} = $id;
		push  @genes, $gene;
		print "gene: $gene, $s, $e\n";
	    }
	    $pre_end = $e;      #need all gene/UTR ends
	}
    }
}
print "Number of genes: $cnt_gene\n";


#parse bam
if($in =~ /\.sam$/){
    open(IN, $in) or die "can't open in sam file: $in\n";
}elsif($in =~ /\.bam$/){
    open(IN, "/hgsc_software/samtools/samtools-1.9/bin/samtools view -h $in|") or die "can't samtools pipe bam: $in\n";
}

open(OUT, ">$out") or die "can't open out: $out\n";
open(OUTSJ, ">$out_sj") or die "can't open splice junction out: $out_sj\n";

my %HOSJ =();   #splicing juntions, index by chr start end
my %HOGR =();   #gene read count

my ($cnt_map_read,$cnt_splice_read_all) = (0,0);
LANE:while(<IN>){
    if(/^\@/){
	#print OUT $_;
    }elsif (/^\S+\s+\d+/){             #match
	chomp;
	my($q,$fl,$t,$ts,$phred,$cig,$matet,$mts,$libsize,$rfa,$rqa,$tags) = split '\t', $_;
	if($cig !~ /\*/){   #cig string can be * if its mate match
	    $cnt_map_read++;
	    my $splice_start = $ts -1;      #5' splicing coord
	    my $splice_end =  '';           #3' splicing
	    #while($cig =~ /(\d+)([MIDNSHP])/g){ #eg: cigar string: 10M1I8M1D56M 1S69M5S
	    while($cig =~ /(\d+)([MIDNPX])/g){ #eg: cigar string: igored 'SH' for read clipping
		my($ln,$type) =($1,$2);

		if($type eq 'D'){  #splicing looks like deletion in reads
		    #print "deletion: $q, $cig: $ln bp\n";
		    if($ln >= $min_splice ){
			if($splice_start >= $leader_start && $splice_start <= $leader_end){  #5' splice coord in leader
			    $splice_end = $splice_start + $ln;
			    $HOSJ{$t}{$splice_start}{$splice_end}++; 
			    $cnt_splice_read_all++;
			    #print OUT $_, "\n";
			    next LANE;
			}
		    }
		}elsif($type eq 'N'){  #N is for splicing for STAR alignment
		    #print "aplicing: $q, $cig: $ln bp\n";
		    if($ln >= $min_splice){
			if($splice_start >= $leader_start && $splice_start <= $leader_end){  #5' splice coord in leader
			    $splice_end = $splice_start + $ln;
			    $HOSJ{$t}{$splice_start}{$splice_end}++;
			    $cnt_splice_read_all++;
			    #print OUT $_, "\n";
			    next LANE;
			}
		    }
		}
		#update target match length for MDNX
		if($type =~ /[MDNX]/){         #note: ISHP type don't change ref coordinates
		    $splice_start += $ln;
		}

	    }
	    
	}
    }
}

close IN;

print "Number of mapped reads: $cnt_map_read\n";
print "Number of splcing reads: $cnt_splice_read_all\n";

#print splice junction, and assign splice read to genes %HOGR
#my $cnt_gene = scalar(@genes);
foreach my $chr (sort keys %HOSJ){
    foreach my $ss (sort {$a <=> $b} keys %{$HOSJ{$chr}} ){
	foreach my $se (sort {$a <=> $b} keys %{$HOSJ{$chr}{$ss}} ){   #3' splice co
	    my $cnt_splice_read = $HOSJ{$chr}{$ss}{$se};
	    #assign splice to gene
	    my $assigned_gene = '.';
	    #my ($pre_gene, $pre_gene_s) = ();    #previous gene, and ATG start coord
	  LANE:for(my $i= 0; $i<$cnt_gene; $i++ ){  #gene are all +, and coordinate sorted ?
	      my $gene = $genes[$i];
	      my $gene_chr = $HOGU{$gene}{'chr'};
	      my $upstream = $HOGU{$gene}{'upstream'};  #upstream cord to end of previous orf
	      my $gene_s  = $HOGU{$gene}{'s'};          #ATG start
	      my $gene_e  = $HOGU{$gene}{'e'};
	      
	      if($gene_chr eq $chr && $se >= $upstream &&  $se <= $gene_s ){  #splice between upstream -> ATG, count this gene
		  $HOGR{$gene} += $cnt_splice_read;
		  $assigned_gene = $gene;
		  last LANE;
	      }elsif($gene_chr eq $chr && $i < ($cnt_gene - 1) && $se >  $gene_s && $se <= $HOGU{$genes[$i+1]}{'s'}){  #assigned to next gene
		  my $next_gene =  $genes[$i+1];
		  my $next_gene_s  = $HOGU{$next_gene}{'s'};
		  my $next_gene_e = $HOGU{$next_gene}{'e'};
		  $HOGR{$next_gene} += $cnt_splice_read;
		  $assigned_gene =  $next_gene;
		  #print "pass start codon: assigned to next gene: $chr,$ss,$se; current_gene: $gene,$gene_s,$gene_e; next_gene: $next_gene,$next_gene_s,$next_gene_e\n";
		  last LANE;
	      }
	      #($pre_gene, $pre_gene_s) = ($gene, $gene_s);
	  }
	    print OUTSJ join("\t", $chr,$ss,$se,$cnt_splice_read,$assigned_gene), "\n";
	}
    }
}


print OUT join("\t","#Gene/Sample", $sample),"\n";
print OUT join("\t","Mapped_reads",$cnt_map_read),"\n";
print OUT join("\t","Junction_reads_with_leader",$cnt_splice_read_all),"\n";

foreach my $gene (@genes){
    my $gene_read_cnt = $HOGR{$gene}?$HOGR{$gene}:'0';
    print OUT join("\t",$gene,$gene_read_cnt), "\n";
}

close OUTSJ;
close OUT;



