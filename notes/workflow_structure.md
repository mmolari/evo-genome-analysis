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

## pileup workflow

For each reference file `{ref_id}` multiple records `{rec_id}` can be present. For each pair of reference and sample `{sample_id}` given in the configuration file the following operations are performed in the workflow:

```mermaid
flowchart TD
	Sample["{sample_id}.fastq.gz"]
	Ref("{ref_id}.fa")
	Map("mapped_reads/{ref_id}/
    {sample_id}.sam,
    {sample_id}.bam,
    {sample_id}.bam.bai")
  Info("mapped_reads/{ref_id}/{sample_id}.info.txt")
  Pileup("pileup/{ref_id}/{rec_id}/{sample_id}/
    allele_counts.npz
    insertions.pkl.gz
    clips.pkl.gz")
  Unmap("pileup/{ref_id}/non_primary/{sample_id}/
    unmapped.csv,
    unmapped.fastq.gz")
  NonPrim("pileup/{ref_id}/non_primary/{sample_id}/
    non_primary.csv")
  Cts("/pileup/{ref_id}/{rec_id}/{coverage,gaps,consensus}.npz")
  Icts("/pileup/{ref_id}/{rec_id}/insertions.npz")
  Ccts("/pileup/{ref_id}/{rec_id}/clips.npz")

	Sample --> |map_reads| Map
	Ref --> Map
  Map --> |map_summary| Info
  Map --> |build_pileup| Pileup
  Map --> |extract_unmapped| Unmap
  Map --> |extract_nonprimary| NonPrim
  Pileup --> |extract_counts| Cts
  Ref --> Cts
  Pileup --> |extract_insertion_counts| Icts
  Ref --> Icts
  Pileup --> |extract_clip_counts| Ccts
  Ref --> Ccts
```

### output file description by rule
- `map_summary`: creates a text file with information on the number and distribution of mapped reads on the various references, for primary, secondary and supplementary alignments.
- `extract_unmapped`: rule to extract unmapped reads. These will contain possible contaminations that do not map to any of the references.
  - `unmapped.csv`: dataframe with one entry per unmapped read. Columns are query name, read length, average quality score and read flag in the sam file.
  - `unmapped.fastq.gz`: fastq file containing unmapped reads.
- `extract_nonprimary`: rule to create a csv file containing info
  - `non_primary.csv`: contains one line per mapping, and the following columns:
    - `qry`/`ref`: reqference and query name
    - `qry_len`: query length. It is the length of the query sequence inferred from the cigar string, including hard clips (it's zero if the sequence is unmapped).
    - `ref_len`: length of the aligned portion on the reference.
    - `n_matches`: total n. of nucleotide matches in the cigar string.
    - `flag`: sam file flag of the mapping.
    - `fwd`/`sec`/`suppl`: whether the mapping is forward / secondary / supplementary
    - `rs`/`re`/`qs`/`qe`: reference/query start and end positions. Differently from the sam file, this is always in the forward frame of reference, so that start and end points can be compared for different mappings.
- `extract_counts`: extract counts of coverage, consensus and gaps. These are (2,L) matrices, where L is the length of the reference sequence and the first index corresponds to fwd/rev reads. They contain integer values indicating counts.
- `extract_insertion_counts`: extract counts of insertion number/length. These are (4,L) matrices - one per sample - contaning per position pair:
  - fwd/rev number of insertions (idx 0-1).
  - total fwd/rev inserted length (idx 2-3).
- `extract_clip_counts`: extract counts of read clip number/length. These are (6,L) matrices - one per sample - contaning per position pair:
  - fwd/rev number of clips (idx 0-1).
  - fwd/rev n. of read ends (idx 2-3).
  - total fwd/rev clipped length (idx 4-5).


# conda environments

The `pileup.yml` environment was created with the command:
```bash
mamba create -n pileup -c bioconda -c conda-forge minimap2 samtools pysam biopython pandas
```

The `plots.yml` environment instead:
```bash
 mamba create -n plots -c conda-forge  numpy scipy matplotlib seaborn pandas plotly
```