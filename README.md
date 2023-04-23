# genome evolution experiment analysis

This pipeline can be used for the genomic analysis of evolution experiments, in which a population is sequenced multiple times over the course of the experiment.

Example run command for local execution:
```bash
snakemake -c1 --profile local --config run_config=test_data/run_config.yml
```

Example run command for cluster execution:
```bash
snakemake --profile cluster --config run_config=test_data/run_config.yml
```
