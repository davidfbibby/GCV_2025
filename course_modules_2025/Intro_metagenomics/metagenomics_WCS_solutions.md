# Solutions

**1.	Write a command to trim adaptors and low quality regions from your data.**
```
trim_galore -q 15 --length 60 --paired ~/metagenomics/data/sample1_1.fq.gz ~/metagenomics/data/sample1_2.fq.gz
```

**2. Write commands to align the reads to the human genome.**
```
bwa index ~/metagenomics/data/human_genome.fasta
bwa mem -t 4 ~/metagenomics/data/human_genome.fasta sample1_1_val_1.fq.gz sample1_2_val_2.fq.gz > ~/metagenomics/sample1.sam
```

**3. What is the output format of the alignment?**

.sam file

**4. How could we use our alignment to get the non-human reads?**

Extract the non aligned reads (in this case to a bam file).

Convert the resulting bam file back to fastq.

**5. Write a command to run kraken2 on your data.**
```
kraken2 --db ~/metagenomics/kraken_db \
--paired ~/metagenomics/sample1_nonhuman_1.fq \
~/metagenomics/sample1_nonhuman_2.fq \
--report ~/metagenomics/sample1_kraken_report.txt \
> ~/metagenomics/sample1_kraken_readclassifications.txt
```

**6. Write a command to run bracken on your kraken2 output.**
```
bracken -d ~/metagenomics/kraken_db \
-i ~/metagenomics/sample1_kraken_report.txt \
-o ~/metagenomics/sample1_bracken.txt \
-t 3
```

**7. What species are present in your sample?**
 
Sample1 contains human mastadenovirus F and cytomegalovirus.

**8. Look at the negative control. What can you conclude about what might be causing the disease in the patient?**
   
The negative control also contains ~5 reads of cytomegalovirus so this is probably a contaminant. Therefore, we would report only the adenovirus clinically.

**9-11.    Extension questions: ask if you need help!**

**12.    Are there any samples that look different to the others? What might these be?**

H20_1 and H20_2 are negative controls.

**13.    Roughly how many raw reads were there for each sample?**

Number of raw reads ranges from around 3 million to 150 million.

**14.    Roughly what proportion of the reads passed quality control? What about filtering? What do these metrics mean?**

During quality control, low quality and complexity and short reads are removed. Filtering happens after QC and is when the human reads are removed. In this dataset, typically 50-90% of reads pass QC and less than 3% pass filtering due to the human human content of the samples.

**15. During which part of QC or filtering were most of the reads lost?**

Most reads were lost during human filtering, followed by low quality filtering.

**16. How can you filter the output to identify species that were present at higher abundances in the sample than in the negative control?**

Filter Z score >= 1.

**17. What is the main viral infection that we might report for this patient? Are there any others?**

Chikungunya virus was found at high levels. Human mastadenovirus C, Rotavirus A and Human alphaherpesvirus 2 are found at lower levels.

**18. Visualise the coverage for the viruses you've found. What do you notice?**

A complete genome with good depth is obtained for Chikungunya virus. For the other viruses, coverage is more patchy. 

**19-21.    Extension questions: ask if you need help!**
