import yaml

config = yaml.safe_load(open("config.yaml"))
FASTA_FILES = [f.split("/")[-1].split(".")[0] for f in config["fasta_files"]]
K_VALUES = config["ks"]
T_VALUES = config["thresholds"]
MODES = ["simplitigs", "unitigs"]  # Ajout des modes

rule all:
    input:
        expand("Benchmark/results_{mode}/{dataset}", mode=MODES, dataset=FASTA_FILES)

rule benchmark:
    input:
        fasta=lambda wildcards: f"test-data/{wildcards.dataset}.fasta"
    output:
        directory("Benchmark/results_{mode}/{dataset}")
    params:
        k_values=K_VALUES,
        t_values=T_VALUES,
        mode=lambda wildcards: wildcards.mode
    shell:
        """
        mkdir -p {output}
        for k in {params.k_values}; do
            for t in {params.t_values}; do
                python src/sequences_to_indexed_spss.py -i {input.fasta} -k $k -t $t --mode {params.mode} -o {output}
            done
        done
        """
