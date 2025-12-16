# workflow/helpers.smk 

def get_best_assembly(wildcards):
    """All downstream analyses use this"""
    return f"output/assemblies/best/{wildcards.sample}.fasta"
