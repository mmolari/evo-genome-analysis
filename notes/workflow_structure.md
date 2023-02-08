# Pipeline setup

## general principles

Separate pileup and plots part.

the **pileup** part is only a function of read and reference (possibly multiple references).

```mermaid
flowchart TD
	A("read_id.fastq.gz")
	B("ref_id.fa")
	C("read_id/ref_id.{sam,bam,bam.bai}")
	D("read_id/ref_id folder")

	A --> |minimap2 + samtools| C
	B --> C
	C --> |python scripts| D
```
The `read_id/ref_id` folder contains
- some summary file on the distribution of primary/secondary reads vs different references (matrix?).
- a set of `record_id` subfolders with the pileup, clips and insertions for reads having primary mapping on the corresponding record, and possibly primary-secondary or primary-supplementary mappings on the same record.

The **plots** part will take these output files to produce summary plots and summary csv files.

## input
```
input/
├── reads
│   ├── sample_1.fastq.gz
│   ├── sample_2.fastq.gz
│   └── ...
└── references
    ├── ref_1.fa
    ├── ref_2.fa
    └── ...
```
**References** must be fasta files (with short record names?)
**Reads** must be fastq files with extension `fastq.gz`

## config file

For the **pileup** part, the config folder contains pairs of reference fasta files and corresponding set of reads to be mapped. Input and output folders can either be specified or have default values.
```yaml
input: "input_fld"
output: "output_fld"
ref_1:
  - "reads_1"
  - "reads_2"
  - ...
ref_2:
  - "reads_1"
  - "reads_2"
  - ...
```

**Plots**: the config file contains set of figures, each one with its own reference and samples set.
```yaml
input: "input_fld"
output: "output_fld"
figure_1:
  reference: "ref_1.fa"
  samples:
    - "reads_1.fastq.gz"
    - "reads_2.fastq.gz"
    - "reads_3.fastq.gz"
    - ...
figure_2:
  reference: "ref_1.fa"
  samples:
    - "reads_1.fastq.gz"
    - "reads_5.fastq.gz"
    - ...
figure_3:
  reference: "ref_2.fa"
  samples:
    - "reads_1.fastq.gz"
    - "reads_6.fastq.gz"
    - ...
```