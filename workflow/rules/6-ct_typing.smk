# Ct Typing Module for ctgap
# Integrates ANI + Mash typing for assembled genomes
# Works with all assembly modes via unified "best" assembly

import pandas as pd
import os
from pathlib import Path

# CT typing configuration
CT_CONFIG = config.get("ct_typing", {})
CT_ENABLED = CT_CONFIG.get("enabled", False)

if CT_ENABLED:
    # Get CT reference genomes from resources
    CT_REF_DIR = RESOURCES / "references" / ORG / "individual"
    CT_REFS = sorted([str(f.resolve()) for f in CT_REF_DIR.glob("*.fa*") 
                      if f.suffix in [".fa", ".fna", ".fasta", ".gz"]])
    
    if not CT_REFS:
        raise ValueError(f"No CT reference genomes found in {CT_REF_DIR}")
    
    # Output directory for CT typing
    CT_OUTDIR = OUTDIR / "ct_typing"
    
    # ========== GLOBAL: Sketch CT references once ==========
    
    rule ct_sketch_refs:
        """Sketch all CT reference genomes once (global step)"""
        input:
            refs = CT_REFS
        output:
            sketch = CT_OUTDIR / "refs" / "ct_refs.msh"
        params:
            k = CT_CONFIG.get("mash", {}).get("k", 31),
            s = CT_CONFIG.get("mash", {}).get("sketch_size", 50000),
            outdir = CT_OUTDIR / "refs"
        threads: 1
        conda: "../envs/ct_typing.yaml"
        log: CT_OUTDIR / "logs" / "ct_sketch_refs.log"
        shell: """
        mkdir -p {params.outdir}
        mash sketch -o {params.outdir}/ct_refs -k {params.k} -s {params.s} {input.refs} 2> {log}
        """
    
    # ========== PER-SAMPLE: Mash distance ==========
    
    rule ct_mash_sketch_query:
        """Sketch individual query genome"""
        input:
            query = get_best_assembly  # ← CHANGED: Use unified best assembly
        output:
            sketch = OUTDIR / "{sample}" / "ct_typing" / "mash" / "{sample}.msh"
        params:
            k = CT_CONFIG.get("mash", {}).get("k", 31),
            s = CT_CONFIG.get("mash", {}).get("sketch_size", 50000),
            prefix = lambda w: OUTDIR / w.sample / "ct_typing" / "mash" / w.sample
        threads: 1
        conda: "../envs/ct_typing.yaml"
        log: OUTDIR / "{sample}" / "log" / "ct_typing.mash_sketch.{sample}.log"
        shell: """
        mkdir -p $(dirname {params.prefix})
        mash sketch -o {params.prefix} -k {params.k} -s {params.s} {input.query} 2> {log}
        """
    
    rule ct_mash_dist:
        """Calculate Mash distances between query and all CT refs"""
        input:
            query_sketch = rules.ct_mash_sketch_query.output.sketch,
            ref_sketch = rules.ct_sketch_refs.output.sketch
        output:
            dist = OUTDIR / "{sample}" / "ct_typing" / "mash" / "{sample}_vs_refs.dist.tsv"
        threads: 1
        conda: "../envs/ct_typing.yaml"
        log: OUTDIR / "{sample}" / "log" / "ct_typing.mash_dist.{sample}.log"
        shell: """
        mash dist {input.ref_sketch} {input.query_sketch} | sort -k3,3g > {output.dist} 2> {log}
        """
    
    rule ct_mash_best:
        """Select best Mash hit for this sample"""
        input:
            dist = rules.ct_mash_dist.output.dist
        output:
            best = OUTDIR / "{sample}" / "ct_typing" / "mash" / "{sample}_best_mash.tsv"
        run:
            df = pd.read_csv(input.dist, sep="\t", header=None,
                           names=["ref", "query", "mash_dist", "pvalue", "shared_hashes"])
            
            # Get best hit (lowest distance)
            best = df.sort_values("mash_dist").iloc[0:1]
            best["sample"] = wildcards.sample
            
            best.to_csv(output.best, sep="\t", index=False)
    
    # ========== PER-SAMPLE: fastANI ==========
    
    rule ct_fastani:
        """Run fastANI between query and all CT references"""
        input:
            query = get_best_assembly,  # ← CHANGED: Use unified best assembly
            refs = CT_REFS
        output:
            ani = OUTDIR / "{sample}" / "ct_typing" / "ani" / "{sample}_vs_refs.ani.tsv",
            ref_list = temp(OUTDIR / "{sample}" / "ct_typing" / "ani" / "refs.txt"),
            query_list = temp(OUTDIR / "{sample}" / "ct_typing" / "ani" / "query.txt")
        params:
            fragLen = CT_CONFIG.get("fastani", {}).get("fragLen", 3000),
            minFraction = CT_CONFIG.get("fastani", {}).get("minFraction", 0.3)
        threads: CT_CONFIG.get("fastani", {}).get("threads", 8)
        conda: "../envs/ct_typing.yaml"
        log: OUTDIR / "{sample}" / "log" / "ct_typing.fastani.{sample}.log"
        shell: """
        mkdir -p $(dirname {output.ani})
        
        # Create file lists
        printf "%s\\n" {input.refs} > {output.ref_list}
        echo "{input.query}" > {output.query_list}
        
        # Run fastANI
        fastANI --rl {output.ref_list} --ql {output.query_list} \
                -o {output.ani} \
                --fragLen {params.fragLen} \
                --minFraction {params.minFraction} \
                --threads {threads} 2> {log} || touch {output.ani}
        """
    
    rule ct_fastani_best:
        """Select best fastANI hit and check for near-ties"""
        input:
            ani = rules.ct_fastani.output.ani,
            mash_all = rules.ct_mash_dist.output.dist
        output:
            best = OUTDIR / "{sample}" / "ct_typing" / "ani" / "{sample}_best_ani.tsv"
        params:
            delta = CT_CONFIG.get("decision", {}).get("delta_ani", 0.001),
            min_cov = CT_CONFIG.get("decision", {}).get("min_ani_cov", 0.95),
            mash_margin = CT_CONFIG.get("decision", {}).get("min_mash_margin", 5e-5),
            need_agree = CT_CONFIG.get("decision", {}).get("require_methods_agree", 2)
        run:
            # Read fastANI results
            if os.path.getsize(input.ani) == 0:
                # No ANI results (assembly too fragmented or different)
                result = pd.DataFrame([{
                    "sample": wildcards.sample,
                    "ref": "NO_ANI_RESULT",
                    "ANI": pd.NA,
                    "frags_over": 0,
                    "frags_total": 0,
                    "cov": 0.0,
                    "alt_ref": "",
                    "alt_ANI": pd.NA,
                    "alt_cov": pd.NA,
                    "flag": "no_ani_result"
                }])
            else:
                df = pd.read_csv(input.ani, sep="\t", header=None,
                               names=["query", "ref", "ANI", "frags_over", "frags_total"])
                
                # Calculate coverage
                df["cov"] = df["frags_over"] / df["frags_total"]
                
                # Sort by ANI, then coverage
                df = df.sort_values(["ANI", "cov"], ascending=[False, False]).reset_index(drop=True)
                
                # Get top hit
                top = df.iloc[0].copy()
                result = pd.DataFrame([{
                    "sample": wildcards.sample,
                    "ref": top["ref"],
                    "ANI": top["ANI"],
                    "frags_over": top["frags_over"],
                    "frags_total": top["frags_total"],
                    "cov": top["cov"],
                    "alt_ref": "",
                    "alt_ANI": pd.NA,
                    "alt_cov": pd.NA,
                    "flag": "clear"
                }])
                
                # Check coverage
                if top["cov"] < params.min_cov:
                    result.loc[0, "flag"] = "low_ani_coverage"
                
                # Check for near-tie with second-best
                if len(df) > 1:
                    alt = df.iloc[1]
                    result.loc[0, "alt_ref"] = alt["ref"]
                    result.loc[0, "alt_ANI"] = alt["ANI"]
                    result.loc[0, "alt_cov"] = alt["cov"]
                    
                    if (top["ANI"] - alt["ANI"]) < params.delta:
                        result.loc[0, "flag"] = "near_tie"
                        
                        # Vote using Mash to break tie
                        mash_df = pd.read_csv(input.mash_all, sep="\t", header=None,
                                            names=["ref", "query", "mash_dist", "pvalue", "shared_hashes"])
                        
                        top_ref = os.path.basename(str(top["ref"]))
                        alt_ref = os.path.basename(str(alt["ref"]))
                        
                        # Extract basename for matching
                        mash_df["ref_base"] = mash_df["ref"].apply(lambda x: os.path.basename(str(x)))
                        
                        m_top = mash_df[mash_df["ref_base"] == top_ref]["mash_dist"].min()
                        m_alt = mash_df[mash_df["ref_base"] == alt_ref]["mash_dist"].min()
                        
                        if pd.notna(m_top) and pd.notna(m_alt):
                            if (m_alt - m_top) > params.mash_margin:
                                result.loc[0, "flag"] = "clear_by_mash_vote"
                            elif (m_top - m_alt) > params.mash_margin:
                                result.loc[0, "ref"] = alt["ref"]
                                result.loc[0, "ANI"] = alt["ANI"]
                                result.loc[0, "cov"] = alt["cov"]
                                result.loc[0, "flag"] = "alt_by_mash_vote"
            
            result.to_csv(output.best, sep="\t", index=False)
    
    # ========== PER-SAMPLE: Final CT assignment ==========
    
    rule ct_assign:
        """Final CT type assignment for individual sample"""
        input:
            mash_best = rules.ct_mash_best.output.best,
            ani_best = rules.ct_fastani_best.output.best
        output:
            assignment = OUTDIR / "{sample}" / "ct_typing" / "{sample}.ct_assignment.csv"
        params:
            high_ani = CT_CONFIG.get("thresholds", {}).get("high_ani", 99.95),
            med_ani = CT_CONFIG.get("thresholds", {}).get("med_ani", 98.0),
            min_cov = CT_CONFIG.get("decision", {}).get("min_ani_cov", 0.95)
        run:
            # Read results
            mash = pd.read_csv(input.mash_best, sep="\t")
            ani = pd.read_csv(input.ani_best, sep="\t")
            
            # Merge
            result = ani.merge(mash.rename(columns={"ref": "ref_mash"}), 
                             on="sample", how="left")
            
            # Extract clean reference names (remove path, keep just filename stem)
            def clean_ref_name(ref_path):
                if pd.isna(ref_path) or ref_path == "NO_ANI_RESULT":
                    return ref_path
                p = Path(str(ref_path))
                # Remove .fa, .fasta, .fna, .gz extensions
                name = p.name
                for ext in [".fa.gz", ".fna.gz", ".fasta.gz", ".fa", ".fna", ".fasta"]:
                    if name.endswith(ext):
                        return name[:-len(ext)]
                return p.stem
            
            result["ref_clean"] = result["ref"].apply(clean_ref_name)
            result["ref_mash_clean"] = result["ref_mash"].apply(clean_ref_name)
            result["alt_ref_clean"] = result["alt_ref"].apply(lambda x: clean_ref_name(x) if x else "")
            
            # Determine confidence
            def get_confidence(row):
                ani_val = row["ANI"]
                cov_val = row["cov"]
                flag = row["flag"]
                
                if pd.isna(ani_val) or flag == "no_ani_result":
                    return "Low/Needs review (no ANI)"
                
                if cov_val < params.min_cov:
                    return "Low/Needs review (low coverage)"
                
                if flag in ["near_tie"]:
                    return "Review (near-tie)"
                
                if ani_val >= params.high_ani:
                    return f"High (ANI≥{params.high_ani}%)"
                elif ani_val >= params.med_ani:
                    return f"Medium ({params.med_ani}–{params.high_ani}%)"
                else:
                    return "Low/Needs review"
            
            result["confidence"] = result.apply(get_confidence, axis=1)
            
            # Choose putative reference (prefer ANI if good coverage)
            result["putative_reference"] = result.apply(
                lambda r: r["ref_clean"] if (pd.notna(r["ANI"]) and r["cov"] >= params.min_cov)
                         else r["ref_mash_clean"],
                axis=1
            )
            
            # Reorder columns for clarity
            cols = [
                "sample", "putative_reference", "confidence",
                "ANI", "cov", "frags_over", "frags_total",
                "ref_clean", "alt_ref_clean", "alt_ANI", "alt_cov", "flag",
                "mash_dist", "pvalue", "shared_hashes", "ref_mash_clean"
            ]
            result = result[[c for c in cols if c in result.columns]]
            
            result.to_csv(output.assignment, index=False)
    
    # ========== COMBINED: Aggregate all samples ==========
    
    rule ct_typing_summary:
        """Aggregate CT typing results for all samples"""
        input:
            assignments = expand(OUTDIR / "{sample}" / "ct_typing" / "{sample}.ct_assignment.csv",
                               sample=SAMPLES)
        output:
            summary = OUTDIR / "reports" / "ct_typing_all_samples.csv",
            status = OUTDIR / "status" / "ct_typing.summary.txt"
        run:
            # Concatenate all sample results
            dfs = []
            for f in input.assignments:
                df = pd.read_csv(f)
                dfs.append(df)
            
            combined = pd.concat(dfs, ignore_index=True)
            
            # Sort by sample name
            combined = combined.sort_values("sample")
            
            # Save
            combined.to_csv(output.summary, index=False)
            
            # Create status file
            Path(output.status).touch()
    
    # ========== Status tracking per sample ==========
    
    rule ct_typing_status:
        """Create status file for each sample's CT typing completion"""
        input:
            assignment = rules.ct_assign.output.assignment
        output:
            status = OUTDIR / "status" / "ct_typing.{sample}.txt"
        shell: """
        touch {output.status}
        """
