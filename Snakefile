import yaml

config = yaml.safe_load(open("config.yaml"))
FASTA_FILES = [f.split("/")[-1].split(".")[0] for f in config["fasta_files"]]
K_VALUES = config["ks"]
T_VALUES = config["thresholds"]

rule all:
    input:
        expand("Benchmark/results/{dataset}/", dataset=FASTA_FILES)

rule benchmark:
    input:
        fasta=lambda wildcards: f"test-data/{wildcards.dataset}.fasta"
    output:
        directory("Benchmark/results/{dataset}/")
    params:
        k_values=K_VALUES,
        t_values=T_VALUES
    shell:
        """
        mkdir -p {output}
        for k in {params.k_values}; do
            for t in {params.t_values}; do
                python src/sequences_to_indexed_spss.py -i {input.fasta} -k $k -t $t -o {output}
            done
        done
        """
