---

<p align="center">
  <img src="./images/protinseq.png"/>
</p>

A transposon sequencing protocol that selects insertions in-frame to expressed genes.

# Workflow 

## 1. Insertion calling 

Once the Tn-Seq data is produced with the modified transposon, these can be normaly processed using your favorite insertion caller to retrieve the genome positions contiguous to the Inverted Repeat (IR) used. 

We recommend the use of [FASTQINS](https://github.com/CRG-CNAG/fastqins). Please follow the previous link for details on the installation of this tools. Keep in mind specific libraries are required by this tool, including standard tools commonly used in high-throughput sequencing analysis:

  [Fastuniq](https://sourceforge.net/projects/fastuniq/) <br /> 
  [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)<br />
  [Samtools](http://www.htslib.org/)<br/>
  [Bedtools](https://bedtools.readthedocs.io/en/latest/)

### Example on how to run FASTQINS for ProTInSeq

Requirements to run an experiment are: 

  -i [fastq files with transposon mapped, if no -i2 is passed, single-end mapping by default] <br />
  -t [IR transposon sequence, expected to be found contiguous genome sequence] <br />
  -g [genome sequence, fasta or genbank format]  <br />
  -o [output directory to locate the results]

As example, a pair of files that you can use to test the pipeline are included in the repository:

```bash
fastqins -i ./test/test_read2.fastq.gz -i2 ./test/test_read1.fastq.gz -t TACGGACTTTATC -g ./test/NC_000912.fna -o test -v -r 0
```

To see additional arguments:
```bash
fastqins --help
```

### Output Information:

The following files are generated as default output:
- \*_fw.qins - read counts of insertions mapping to forward strand \[[example](./test/output_test/test_read2_fw.qins)\]
- \*_rv.qins - read counts of insertions mapping to reverse strand \[[example](./test/output_test/test_read2_rv.qins)\]
- \*.qins - read counts of insertions mapping to both strands \[[example](./test/output_test/test_read2.qins)\]
- \*.bam - file generated with the aligned reads
- \*.log - log file with general features of the process run \[[example](./test/output_test/test_read2.log)\]

## 2. Analysis of the insertion profiles

For control libraries, no difference is expected. Please refer to our previous publication [FASTQINS and ANUBIS: two bioinformatic tools to explore facts and artifacts in transposon sequencing and essentiality studies](https://academic.oup.com/nar/article/48/17/e102/5894413) to follow the best practices in analyzing this data. 

We recommend [ANUBIS](https://github.com/CRG-CNAG/anubis) for the analysis of this type of data. However, simple data analysis can be performed to extract a relation of genomic annotations and relevant metrics by frame. We have included a set of useful functions in the file [protinseq.py](./protinseq.py). 

A small demonstration on how to apply them is included in [demonstration](./protinseq_analysis.ipynb). 

Please ensure you have the most recent installation of the packages in the [requirements](./requirements.txt) by:
```bash
pip install -r requirements.txt
```

### Output examples:

By running the previous notebook you can obtain:

- Comprehensive omics information for *M. pneumoniae*.
- Sample exploration.
- Table with metrics by genomic annotation including genes, smORFs and intergenic regions (used as control)
- Basic plotting of loci of interest coloring insertions by frame.
- Metagene exploration.

# System requirements

No special system requirements are required to run this pipelines. The presented analysis has been run in a Linux operative system and tested in MacOS and Windows running WSL. We expect a Python 3.6 or higher version to run the processes. 

# Contact

This project has been fully developed at [Centre for Genomic Regulation](http://www.crg.eu/) at the group of [Design of Biological Systems](http://www.crg.eu/en/luis_serrano).

If you experience any problem at any step involving the program, you can use the 'Issues' page of this repository or contact:

[Miravet-Verde, Samuel](mailto:smiravet@ethz.ch)         
[Serrano, Luis](mailto:luis.serrano@crg.eu)

# Citation

If you use the tools and workflow presented in this repository, please cite:

- [FASTQINS and ANUBIS: two bioinformatic tools to explore facts and artifacts in transposon sequencing and essentiality studies](https://academic.oup.com/nar/article/48/17/e102/5894413)
- ProTInSeq: transposon insertion tracking by ultra-deep DNA sequencing to identify translated large and small ORFs

# License

ProTInSeq is under a common GNU GENERAL PUBLIC LICENSE. Plese, check [LICENSE](./LICENSE) for further information.

###### [2024] - Centre de Regulació Genòmica (CRG) - All Rights Reserved*

