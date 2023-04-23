# input data preparation

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

The name of records in reference files must be short (less than 20 characters) and *not contain any spaces* since it will be used to generate names of the subfiles.

## pileup configuration file

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

this file is passed as input in the pipeline as `run_config` parameter using the `--config` flag:
```bash
snakemake --config run_config=myconfig.yml
```

Moreover the general workflow `config.yaml` file more options can be controlled:
```yaml
pileup:
    qual_min: 15
    clip_minL: 100
```
- `qual_min` controls the minimum quality of reads to be kept in the pileup.
- `clip_minL` controls the minimum length of clips to be saved in the clipped reads dictionary.