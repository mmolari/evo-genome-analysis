# genome evolution experiment analysis

This pipeline can be used for the genomic analysis of evolution experiments, in which a population is sequenced multiple times over the course of the experiment.

Example run command for local execution:
```bash
snakemake --profile local --configfile test_data/run_config.yml
```

Example run command for cluster execution:
```bash
snakemake --profile cluster --configfile test_data/run_config.yml
```

A description of the produced plots is provided in this [note](notes/plot_description.md).

> [!IMPORTANT]  
> The pipeline was developed for snakemake v7. Version 8 might introduce breaking changes for which the pipeline was not tested.