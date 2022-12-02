# VCFtoNPZ

## Purpose
The purpose of this code is to pull .vcf data into a custom numpy format. 

VCF files are compressed text files, often with an index (accompanying .tbi files). VCF files can be extremely large when working with WGS data from many individuals. They contain a huge amount of data including quality scores, allele frequencies, etc, and they take a lot of time to read in. However, most of my work requires only the variant calls for each individual at each position. Since most individuals are homozygous reference most of the time, this data can be stored in a sparse numpy matrix, greatly reducing the memory required to store the data, as well as the compute time needed to parse the VCF file. This is particularly valuable when you have several algorithms you'd like to run on a large WGS dataset.


## Input and output
This code starts with VCF files (with accompanying .tbi files), split by chromosome. It then produces numpy sparse matrix files containing the variant call for every individual at every position. The entries of this sparse matrix are

- -1 represents ./. variant call
- 0 represents 0/0 variant call
- 1 represents 0/1 variant call
- 2 represents 1/1 variant call

For ease of manipulation, the code gives the option to separate each chromosome into batches by genomic position. For example, files can be broken up such that each file represents a 10Mbp chunk of the chromosome.

The code populates a directory (which we call `data_dir` in this documentation) with the following structure.

```
data_dir
- genotypes
- - info.json
- - samples.json
- - chr.1.0.gen.npz
- - chr.1.0.gen.coordinates.npy
- - chr.1.0.gen.variants.txt.gz
- - chr.2.0.gen.npz
- - chr.2.0.gen.coordinates.npy
- - chr.2.0.gen.variants.txt.gz
...
```

The `info.json` file contains metadata including the reference assembly (GRch37 or GRch38) used to produce the variant calls, the vcf_directory used to generate the data, and if relevant the batch_size used to generate the data.

The `samples.json` file contains an ordered list of the samples contained in the dataset.

The `chr.[chrom].[batch_num].gen.npz` files contains variant calls with rows representing individuals and columns representing genomic positions. This data is stored in sparse numpy format. The number of rows in each `chr.[chrom].[batch_num].gen.npz` file corresponds to the number of entries in the `samples.json` ordered list.

The `chr.[chrom].[batch_num].gen.coordinates.npy` files contain basic information for each position stored in numpy format. Each row represents a genomic position, with the first column indicating the chromosome, the second column indicating the position (pulled from the POS field of the VCF), the third column indicating whether the variant is a biallelic SNP (1 if yes, 0 if no) and the fourth column indicating whether the variant contains PASS in the FILTER field of the VCF (1 if yes, 0 if no). The number of rows in each `chr.[chrom].[batch_num].gen.coordinates.npy` file corresponds to the number of columns in each `chr.[chrom].[batch_num].gen.npz` file.

The `chr.[chrom].[batch_num].gen.variants.txt.gz` contains the first 8 columns of the vcf file for each position. These include CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO in compressed text format. Once uncompressed, the number of rows in each `chr.[chrom].[batch_num].gen.variants.txt.gz` file corresponds to the number of columns in each `chr.[chrom].[batch_num].gen.npz` file.


## Instructions for running code
1. Create an output directory. Then create a subdirectory within that directory named genotypes. 

2. Run `python preprocessing/pull_gen_data.py [vcf_file] [assembly] [output_dir] [chrom]`

If you have a large WGS dataset, it may be more convenient to break the data into batches using the `--batch_size` and `--batch_num` flags

3. If your VCF files don't have filters applied (for example no variant is PASS) or you'd like to apply a different type of filter for downstream analysis, use 

`python preprocessing/pull_pass.py [data_dir (output directory from above command)]`

which has options 
`--pass_from_gen [ped_file]` (passes variants on autosomes with <10% missing, passes variants on Y with <20% missing in males, passes variants on X with <10% missing in females and <20% missing in males)
`--pass_all` (passes all variants)
`--pass_from_qual [cutoff]` (passes variants whose QUAL score is better than cutoff)
