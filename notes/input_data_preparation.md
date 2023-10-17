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

## configuration file

The `run_config` entry in the config file is used to instruct the pipeline on where the input data are located, where the output data should be saved, and which reference should be used for each sample.
This must have the following structure:
```yaml
run_config:
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
...
```

Thre config file is then passed as input to the pipeline using the `--configfile` flag:
```bash
snakemake --configfile myconfig.yml
```

The same config file also contains other options that control various parameters of the pipeline. These are described in the [plot description](plot_description.md) file.