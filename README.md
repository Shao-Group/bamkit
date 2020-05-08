# Overview 
bamkit is a toolkit designed to analyze and process `bam/sam` files.
More functions will be continuously added.

# Installation
To install bamkit from the source code, you need to first download/compile 
htslib, setup the corresponding environmental variables,
and then compile the source code of bamkit.


## Install htslib
Download htslib [(license)](https://github.com/samtools/htslib/blob/develop/LICENSE)
from (http://www.htslib.org/) with version 1.5 or higher.
Choose a directory for installation and set the environment variable `HTSLIB` for that.
For example, for Unix platforms, do the following:
```
export HTSLIB="/directory/to/your/htslib/install"
```
Use the following to build `libhts.a`:
```
autoheader
autoconf
./configure --disable-bz2 --disable-lzma --disable-gcs --disable-s3 --enable-libcurl=no --prefix=$HTSLIB
make
make install
```

## Compile bamkit
The compilation of `bamkit` requires `automake` and `autoconf` packages.
If they have not been installed, on linux platform, do the following:
```
sudo apt-get install autoconf
sudo apt-get install automake
```

bamkit might also requires other libraries, such as `libz`, depending on
your system. Install them if you encounter errors when compiling.

After that run the script `build.sh`, which will generate the executable file `src/src/bamkit`.


# Usage

The current version of `bamkit` supports the following functionalities:
```
./bamkit ts2XS <input.bam> <output.bam>
```
The `ts` tag (used in `minimap2` aligner) will be transformed into `XS` tag (used in `STAR`, `HISAT` alingers)
	and the resuting alignments will be written to `output.bam`.

```
./bamkit count <input.bam>
```
A statistic will be returned about the `input.bam`, including
the number of reads and basepairs aligned, etc.

```
./bamkit strand <input.bam>
```
A report about the strandness of the `input.bam` will be written to standard output.

```
./bamkit alignPairEval <input.bam> <groundTruth.bam>
```
This command is about evaluation of aligners(like [STAR](https://github.com/alexdobin/STAR)). `input.bam` is the result of aligner and `groundTruth.bam` is the ground truth based on output of [flux simulator](http://confluence.sammeth.net/display/SIM/Home). Evaluation results of the aligner will be written to standard output.


```
./bamkit bridgeEval <input.coral.bam> <input.aligner.bam> <groundTruth.bam> <reference.gtf>
```
This command is about evaluation of aligners and tools that bridge paired reads(like [coral](https://github.com/Shao-Group/coral)). `input.coral.bam` is the result of coral, `input.aligner.bam` is the result of aligner and `groundTruth.bam` is the ground truth based on output of flux simulator. `reference.gtf` is the annotation used to simulate reads. Evaluation results of the aligner and coral will be written to standard output.
