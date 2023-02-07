# input preparation

The input folder must have the following structure:

```
input_fld/
├── reads
│   ├── sample_1.fastq.gz
│   ├── sample_2.fastq.gz
│   └── ...
└── references
    ├── ref_1.fa
    ├── ref_2.fa
    └── ...
```

The name of the reads and reference files is arbitrary, but the extensions must be `.fastq.gz` and `.fa` respectively.

The name of records in reference files must be short (less than 20 characters).

The configuration file for the pileup must have the structure:
```yaml
input: "input_fld"
output: "output_fld"
pileups:
    ref_1:
    - "sample_1"
    - "sample_2"
    - ...
    ref_2:
    - "sample_1"
    - "sample_2"
    - ...
```

this part of the pipeline can be run with:
```bash
snakemake --config run_config=myconfig.yml
```