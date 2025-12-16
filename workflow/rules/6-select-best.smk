# workflow/rules/6-select-best.smk

mode = config.get("mode", "denovo")

if mode == "auto":
    # Run both assemblies, score them, pick winner
    rule select_best_assembly:
        input:
            denovo = OUTDIR / "{sample}" / "denovo" / "denovo_{sample}.final.fasta",  # ← FIXED
            refguided = OUTDIR / "{sample}" / "ref-denovo" / "ref-denovo_{sample}.final.fasta",  # ← FIXED
            denovo_quast = OUTDIR / "{sample}" / "denovo" / "assembly_statistics" / "report.tsv",  # ← FIXED
            refguided_quast = OUTDIR / "{sample}" / "ref-denovo" / "quast" / "report.tsv"  # ← FIXED
        output:
            best = OUTDIR / "{sample}" / "best" / "assembly.fasta",
            report = OUTDIR / "reports" / "assembly_selection" / "{sample}.txt",
            metadata = OUTDIR / "reports" / "assembly_selection" / "{sample}.method.txt"
        log:
            OUTDIR / "logs" / "select_best" / "{sample}.log"
        script:
            "../scripts/select_best_assembly.py"

elif mode == "denovo":
    # Just use denovo assembly
    rule copy_denovo_as_best:
        input:
            assembly = OUTDIR / "{sample}" / "denovo" / "denovo_{sample}.final.fasta"  # ← FIXED
        output:
            best = OUTDIR / "{sample}" / "best" / "assembly.fasta",
            metadata = OUTDIR / "reports" / "assembly_selection" / "{sample}.method.txt"
        shell:
            """
            mkdir -p $(dirname {output.best})
            mkdir -p $(dirname {output.metadata})
            cp {input.assembly} {output.best}
            echo "denovo" > {output.metadata}
            """

elif mode == "reference-denovo":
    # Just use reference-guided assembly
    rule copy_refguided_as_best:
        input:
            assembly = OUTDIR / "{sample}" / "ref-denovo" / "ref-denovo_{sample}.final.fasta"  # ← FIXED
        output:
            best = OUTDIR / "{sample}" / "best" / "assembly.fasta",
            metadata = OUTDIR / "reports" / "assembly_selection" / "{sample}.method.txt"
        shell:
            """
            mkdir -p $(dirname {output.best})
            mkdir -p $(dirname {output.metadata})
            cp {input.assembly} {output.best}
            echo "reference-denovo" > {output.metadata}
            """

else:
    raise ValueError(f"Invalid mode '{mode}'. Must be: auto, denovo, or reference-denovo")


# Summary rule - runs after all samples processed
rule assembly_method_summary:
    input:
        methods = expand(OUTDIR / "reports" / "assembly_selection" / "{sample}.method.txt",
                        sample=SAMPLES)
    output:
        summary = OUTDIR / "reports" / "assembly_methods_summary.txt"
    run:
        from collections import Counter
        methods = []
        for f in input.methods:
            with open(f) as fh:
                methods.append(fh.read().strip())

        counts = Counter(methods)

        with open(output.summary, 'w') as out:
            out.write("Assembly Method Selection Summary\n")
            out.write("=================================\n\n")

            if config["mode"] == "auto":
                out.write("Mode: auto (automatic selection)\n\n")
                out.write(f"Total samples: {len(methods)}\n")
                out.write(f"De novo selected: {counts.get('denovo', 0)} samples\n")
                out.write(f"Reference-guided selected: {counts.get('reference-denovo', 0)} samples\n")
            else:
                out.write(f"Mode: {config['mode']} (single method)\n")
                out.write(f"All {len(methods)} samples used {config['mode']}\n")
