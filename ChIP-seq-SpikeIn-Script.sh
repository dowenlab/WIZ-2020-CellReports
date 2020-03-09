
## ChIP-seq Spike In Normalization Script - Justice et al 2020; Cell Reports
## Author: Megan Justice
## Contact: jilldowen@unc.edu

module load bowtie
module load bedtools
module load samtools
module load fastqc

#A custom genome containing both mm10 and hg38 coordinates was created using bowtie-build. Mouse chromosomes are denoted by an "M" before the chromosome information (i.e. Mchr1). Genome_DIR is the directory to the merged genome files. 

Genome_DIR=/path/to/MergedGenome
chromSizes='/path/to/mm10.chrom.sizes'

## NOTE: For replicates, I merged the fastq files by doing cat file1.fastq file2.fastq before beginning this script. 
DIR=/path/to/fastqFiles


for f in ${DIR}

      do

## Perform quality control and align the fastq files.
      fastqc ${f}
      (bowtie -v 2 -p 24 -S -m 1 --best --strata ${Genome_DIR} ${f} > ${f}.sam) 2>${f}.log

## The sam file must be converted to a sorted bam file.
      samtools view -b -S -F 4 -@ 24 -o ${f}.bam ${f}.sam
      samtools sort -@ 24 -o ${f}-sort.bam ${f}.bam

## Next, remove PCR duplicates.
      samtools sort -n -@ 24 ${f}-sort.bam -o ${f}-sortedbyname.bam
      samtools fixmate -m -@ 24 ${f}-sortedbyname.bam ${f}-fixmate.bam
      samtools sort -@ 24 ${f}-fixmate.bam -o ${f}-fixmate_sorted.bam
      samtools markdup -r -s -@ 24 ${f}-fixmate_sorted.bam ${f}-markdup.bam
      samtools sort -@ 24 -o ${f}-markdup-sort.bam ${f}-markdup.bam
      samtools index ${f}-markdup-sort.bam ${f}-markdup-sort.bam.bai


## Separate the mouse and human reads. Count them and calculate a normalization factor based on Orlando et al 2014. Formula is 1 / # human reads in millions.
      samtools idxstats ${f}-markdup-sort.bam | grep -v Mchr > ${f}-humanalign.txt
      samtools idxstats ${f}-markdup-sort.bam | grep Mchr > ${f}-mousealign.txt
      awk '
              {s+=$3}END{print s}
          ' ${f}-humanalign.txt > ${f}-humancount.txt
      awk '
              {s+=$3}END{print s}
          ' ${f}-mousealign.txt > ${f}-mousecount.txt
	normfactor=$( awk '{h+=$1}END{if(h<50000){print 1} else {print (1 / (h/1000000))}}' ${f}-humancount.txt )

echo $normfactor > ${f}-NormFactor.txt 


## Isolate the sample reads from the spike in reads. Convert to bed format.
      samtools view -h -@ 24 ${f}-markdup-sort.bam | grep Mchr > ${f}-mousereads.bam
      sed 's/Mchr/chr/g' ${f}-mousereads.bam > ${f}-mousereads-fix.bam
      samtools sort ${f}-mousereads-fix.bam > ${f}-mousereads-fix-sort.bam
      bedtools bamtobed -i ${f}-mousereads-fix-sort.bam > ${f}.bed


## Expand the bed coordinates by 200 bp, ensuring the new coordinates are not off the ends of the chromosome.. This helps smooth the data.
      awk '{

       if ($6 == "-")
               print $1, $2-200, $3, $4, $5, $6;

       else
               print $1, $2, $3+200, $4, $5, $6;

       }' ${f}.bed > ${f}-sort-extend.bed

       awk '{
       if ($2 < 0)
               print $1, 0, $3, $4, $5, $6;
       else
               print $1, $2, $3, $4, $5, $6;

       }' ${f}-sort-extend.bed > ${f}-fixstart.bed

       awk '{OFS="\t"}{
       if (($1 == "chr1") &&
      ($3 > 195471971))

               print $1, $2, 195471971, $4, $5, $6;

     else if (($1 == "chr2") &&
               ($3 > 182113224))

                       print $1, $2, 182113224, $4, $5, $6;

     else if (($1 == "chr3") &&
               ($3 > 160039680))

                       print $1, $2, 160039680, $4, $5, $6;

     else if (($1 == "chr4") &&
               ($3 > 156508116))

                       print $1, $2, 156508116, $4, $5, $6;

     else if (($1 == "chr5") &&
               ($3 > 151834684))

                       print $1, $2, 151834684, $4, $5, $6;

     else if (($1 == "chr6") &&
               ($3 > 149736546))

                       print $1, $2, 149736546, $4, $5, $6;

     else if (($1 == "chr7") &&
               ($3 > 145441459))

                       print $1, $2, 145441459, $4, $5, $6;

     else if (($1 == "chr8") &&
               ($3 > 129401213))

                       print $1, $2, 129401213, $4, $5, $6;
     else if (($1 == "chr9") &&
               ($3 > 124595110))

                       print $1, $2, 124595110, $4, $5, $6;

      else if (($1 == "chr10") &&
               ($3 > 130694993))

                       print $1, $2, 130694993, $4, $5, $6;

      else if (($1 == "chr11") &&
               ($3 > 122082543))

                       print $1, $2, 122082543, $4, $5, $6;

      else if (($1 == "chr12") &&
               ($3 > 120129022))

                       print $1, $2, 120129022, $4, $5, $6;

      else if (($1 == "chr13") &&
               ($3 > 120421639))

                       print $1, $2, 120421639, $4, $5, $6;

      else if (($1 == "chr14") &&
               ($3 > 124902244))

                       print $1, $2, 124902244, $4, $5, $6;

      else if (($1 == "chr15") &&
               ($3 > 104043685))

                       print $1, $2, 104043685, $4, $5, $6;

      else if (($1 == "chr16") &&
               ($3 > 98207768))

                       print $1, $2, 98207768, $4, $5, $6;

      else if (($1 == "chr17") &&
      	($3 > 94987271))

                       print $1, $2, 94987271, $4, $5, $6;

      else if (($1 == "chr18") &&
               ($3 > 90702639))

                       print $1, $2, 90702639, $4, $5, $6;

      else if (($1 == "chr19") &&
               ($3 > 61431566))

                       print $1, $2, 61431566, $4, $5, $6;

      else if (($1 == "chrX") &&
               ($3 > 171031299))

                       print $1, $2, 171031299, $4, $5, $6;

      else if (($1 == "chrY") &&
               ($3 > 91744698))

                       print $1, $2, 91744698, $4, $5, $6;

      else if (($1 == "chrM") &&
               ($3 > 16299))

                       print $1, $2, 16299, $4, $5, $6;



      else if (($1 == "chr5_JH584299_random") &&
               ($3 > 953012))

                       print $1, $2, 953012, $4, $5, $6;

      else if (($1 == "chrX_GL456233_random") &&
              ($3 > 336933))

                       print $1, $2, 336933, $4, $5, $6;
     else if (($1 == "chrY_JH584301_random") &&
               ($3 > 259875))

                       print $1, $2, 259875, $4, $5, $6;

      else if (($1 == "chr1_GL456211_random") &&
               ($3 > 241735))

                       print $1, $2, 241735, $4, $5, $6;

      else if (($1 == "chr4_GL456350_random") &&
               ($3 > 227966))

                       print $1, $2, 227966, $4, $5, $6;

       else if (($1 == "chr4_JH584293_random") &&
               ($3 > 207968))

                       print $1, $2, 207968, $4, $5, $6;

       else if (($1 == "chr1_GL456221_random") &&
               ($3 > 206961))

                       print $1, $2, 206961, $4, $5, $6;

       else if (($1 == "chr5_JH584297_random") &&
               ($3 > 205776))

                       print $1, $2, 205776, $4, $5, $6;

       else if (($1 == "chr5_JH584296_random") &&
               ($3 > 199368))

                       print $1, $2, 199368, $4, $5, $6;

       else if (($1 == "chr5_GL456354_random") &&
               ($3 > 195993))

                       print $1, $2, 195993, $4, $5, $6;

       else if (($1 == "chr4_JH584294_random") &&
               ($3 > 191905))
      		print $1, $2, 191905, $4, $5, $6;

       else if (($1 == "chr5_JH584298_random") &&
               ($3 > 184189))

                       print $1, $2, 184189, $4, $5, $6;

       else if (($1 == "chrY_JH584300_random") &&
               ($3 > 182347))

                       print $1, $2, 182347, $4, $5, $6;

       else if (($1 == "chr7_GL456219_random") &&
               ($3 > 175968))

                       print $1, $2, 175968, $4, $5, $6;

       else if (($1 == "chr1_GL456210_random") &&
               ($3 > 169725))

                       print $1, $2, 169725, $4, $5, $6;

       else if (($1 == "chrY_JH584303_random") &&
               ($3 > 158099))

                       print $1, $2, 158099, $4, $5, $6;

       else if (($1 == "chrY_JH584302_random") &&
               ($3 > 155838))

                       print $1, $2, 155838, $4, $5, $6;

       else if (($1 == "chr1_GL456212_random") &&
               ($3 > 153618))

                       print $1, $2, 153618, $4, $5, $6;

       else if (($1 == "chrUn_JH584304") &&
               ($3 > 114452))

                       print $1, $2, 114452, $4, $5, $6;
      else if (($1 == "chrUn_GL456379") &&
               ($3 > 72385))

                       print $1, $2, 72385, $4, $5, $6;

       else if (($1 == "chr4_GL456216_random") &&
               ($3 > 66673))

                       print $1, $2, 66673, $4, $5, $6;

       else if (($1 == "chrUn_GL456393") &&
               ($3 > 55711))

                       print $1, $2, 55711, $4, $5, $6;

       else if (($1 == "chrUn_GL456366") &&
               ($3 > 47073))

                       print $1, $2, 47073, $4, $5, $6;

       else if (($1 == "chrUn_GL456367") &&
               ($3 > 42057))

                       print $1, $2, 42057, $4, $5, $6;

       else if (($1 == "chrUn_GL456239") &&
               ($3 > 40056))

                       print $1, $2, 40056, $4, $5, $6;

       else if (($1 == "chr1_GL456213_random") &&
               ($3 > 39340))

                       print $1, $2, 39340, $4, $5, $6;

       else if (($1 == "chrUn_GL456383") &&
               ($3 > 38659))

                       print $1, $2, 38659, $4, $5, $6;

       else if (($1 == "chrUn_GL456385") &&
               ($3 > 35240))
      		print $1, $2, 35240, $4, $5, $6;

       else if (($1 == "chrUn_GL456360") &&
               ($3 > 31704))

                       print $1, $2, 31704, $4, $5, $6;

       else if (($1 == "chrUn_GL456378") &&
               ($3 > 31602))

                       print $1, $2, 31602, $4, $5, $6;

       else if (($1 == "chrUn_GL456389") &&
               ($3 > 28772))

                       print $1, $2, 28772, $4, $5, $6;

       else if (($1 == "chrUn_GL456372") &&
               ($3 > 28664))

                       print $1, $2, 28664, $4, $5, $6;

       else if (($1 == "chrUn_GL456370") &&
               ($3 > 26764))

                       print $1, $2, 26764, $4, $5, $6;

       else if (($1 == "chrUn_GL456381") &&
               ($3 > 25871))

                       print $1, $2, 25871, $4, $5, $6;

       else if (($1 == "chrUn_GL456387") &&
               ($3 > 24685))

                       print $1, $2, 24685, $4, $5, $6;

       else if (($1 == "chrUn_GL456390") &&
               ($3 > 24668))

                       print $1, $2, 24668, $4, $5, $6;
      else if (($1 == "chrUn_GL456394") &&
               ($3 > 24323))

                       print $1, $2, 24323, $4, $5, $6;

       else if (($1 == "chrUn_GL456392") &&
               ($3 > 23629))

                       print $1, $2, 23629, $4, $5, $6;

       else if (($1 == "chrUn_GL456382") &&
               ($3 > 23158))

                       print $1, $2, 23158, $4, $5, $6;

       else if (($1 == "chrUn_GL456359") &&
               ($3 > 22974))

                       print $1, $2, 22974, $4, $5, $6;

       else if (($1 == "chrUn_GL456396") &&
               ($3 > 21240))

                       print $1, $2, 21240, $4, $5, $6;

       else if (($1 == "chrUn_GL456368") &&
               ($3 > 20208))

                       print $1, $2, 20208, $4, $5, $6;

       else if (($1 == "chr4_JH584292_random") &&
               ($3 > 14945))

                       print $1, $2, 14945, $4, $5, $6;

       else if (($1 == "chr4_JH584295_random") &&
               ($3 > 1976))

                       print $1, $2, 1976, $4, $5, $6;

       else
               print $1, $2, $3, $4, $5, $6;
      }' ${f}-fixstart.bed > ${f}-fixall.bed


## Call peaks from the expanded bed file. Then, remove peaks that overlap repetitive elements in the genome. RepeatSequence.bed refers to the RepeatMasker track on UCSC's genome browser. This track is available for download via the UCSC Table Browser.
module load macs
module load bedtools       
       macs2 callpeak -g mm -q 0.01 -f BED -n ${f} -t ${f}-fixall.bed
       awk 'BEGIN{OFS="\t"}{print $1, $2-50, $3+50, $4, $5
       }' ${f}_summits.bed > ${f}-summits-exp50.bed
       bedtools intersect -a ${f}-summits-exp50.bed -b RepeatSequence.bed -v -wa > ${f}-summits-exp50-noreps.bed
       awk '{print $1, $2+50, $3-50, $4, $5
       }' ${f}-summits-exp50-noreps.bed > ${f}-summits-exp50-noreps-shrink.bed


## Make a coverage track for the data, scaling according to the normalization factor calculated above.
       bedtools genomecov -bga -g ${chromSizes}  -scale ${normfactor} -i ${f}-fixall.bed > ${f}-fixall.bedgraph


       sort -k1,1 -k2,2n ${f}-fixall.bedgraph > ${f}-sort.bedgraph

module load ucsctools
	bedGraphToBigWig ${f}-sort.bedgraph ${chromSizes} ${f}.bw
done

