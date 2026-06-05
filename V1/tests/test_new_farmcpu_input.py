#!/usr/bin/env python3
import csv
import math
import subprocess
import sys
import tempfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]


def test_rmvp_signals_are_converted_to_farmcpu_peak() -> None:
    script = REPO_ROOT / "V1" / "expand" / "Select_FarmCPU_Peak.py"

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        signals_dir = root / "data_flowering"
        signals_dir.mkdir()
        with (signals_dir / "HNHZ.FarmCPU_signals.csv").open("w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["SNP", "CHROM", "POS", "REF", "ALT", "MAF", "Effect", "SE", "HNHZ.FarmCPU"])
            writer.writerow(["snp1", "1", "1000", "A", "G", "0.1", "0.2", "0.03", "1e-10"])
            writer.writerow(["snp2", "1", "1200", "C", "T", "0.1", "0.4", "0.04", "1e-12"])
            writer.writerow(["snp3", "1", "3005000", "G", "A", "0.2", "-0.5", "0.05", "1e-9"])

        subprocess.run([sys.executable, str(script), "HNHZ"], cwd=root, check=True)

        rows = (root / "HNHZ.FarmCPU.peak").read_text(encoding="utf-8").splitlines()

    assert rows[0] == "Chromosome\tPosition\t-log10(Pvalue)\tParent_geno\t1_genotype\t2_genotype\tEffect\tSE"
    assert rows[1].split("\t") == [
        "1",
        "1200",
        f"{-math.log10(1e-12):.2f}",
        "*",
        "C",
        "T",
        "0.4",
        "0.04",
    ]
    assert rows[2].split("\t") == [
        "1",
        "3005000",
        f"{-math.log10(1e-9):.2f}",
        "*",
        "G",
        "A",
        "-0.5",
        "0.05",
    ]
    assert len(rows) == 3


def test_data_tgw_farmcpu_csv_is_filtered_to_default_significant_sites() -> None:
    script = REPO_ROOT / "V1" / "expand" / "Select_FarmCPU_Peak.py"

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        tgw_dir = root / "data_TGW"
        tgw_dir.mkdir()
        with (tgw_dir / "all_TGW_quality.FarmCPU.csv").open("w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["SNP", "CHROM", "POS", "REF", "ALT", "MAF", "Effect", "SE", "all_TGW_quality.FarmCPU"])
            writer.writerow(["snp_not_significant", "1", "1000", "A", "G", "0.1", "0.2", "0.03", "1e-4"])
            writer.writerow(["snp_at_default", "1", "1200", "C", "T", "0.1", "0.4", "0.04", "1e-5"])
            writer.writerow(["snp_best_peak", "1", "1500", "G", "A", "0.2", "-0.5", "0.05", "1e-6"])
            writer.writerow(["snp_second_peak", "1", "3005000", "T", "C", "0.2", "-0.7", "0.06", "1e-7"])

        subprocess.run([sys.executable, str(script), "all_TGW_quality"], cwd=root, check=True)

        significant = (root / "all_TGW_quality.FarmCPU.significant.tsv").read_text(encoding="utf-8").splitlines()
        peaks = (root / "all_TGW_quality.FarmCPU.peak").read_text(encoding="utf-8").splitlines()

    assert significant[0] == "SNP\tChromosome\tPosition\tPvalue\t-log10(Pvalue)\tREF\tALT\tEffect\tSE"
    assert [row.split("\t")[0] for row in significant[1:]] == [
        "snp_at_default",
        "snp_best_peak",
        "snp_second_peak",
    ]
    assert peaks[1].split("\t")[:3] == ["1", "1500", f"{-math.log10(1e-6):.2f}"]
    assert peaks[2].split("\t")[:3] == ["1", "3005000", f"{-math.log10(1e-7):.2f}"]
    assert len(peaks) == 3


def test_mlm_peak_uses_named_pvalue_column_when_maf_is_present() -> None:
    script = REPO_ROOT / "V1" / "expand" / "Select_MLM_Peak.py"

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        mlm_dir = root / "data" / "farmCPU" / "HNHZ.MLM"
        mlm_dir.mkdir(parents=True)
        with (mlm_dir / "HNHZ.MLM.csv").open("w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["SNP", "CHROM", "POS", "REF", "ALT", "MAF", "Effect", "SE", "HNHZ.MLM"])
            writer.writerow(["snp1", "2", "1000", "A", "G", "0.1", "0.2", "0.03", "1e-9"])
            writer.writerow(["snp2", "2", "1500", "C", "T", "0.1", "0.4", "0.04", "1e-11"])

        subprocess.run([sys.executable, str(script), "HNHZ"], cwd=root, check=True)

        rows = (root / "HNHZ.MLM.peak").read_text(encoding="utf-8").splitlines()

    assert len(rows) == 1
    parts = rows[0].split("\t")
    assert parts[:3] == ["snp2", "2", "1500"]
    assert math.isclose(float(parts[3]), -math.log10(1e-11))


def test_mlm_signals_from_data_flowering_are_supported() -> None:
    script = REPO_ROOT / "V1" / "expand" / "Select_MLM_Peak.py"

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        mlm_dir = root / "data_flowering"
        mlm_dir.mkdir()
        with (mlm_dir / "jap_Ganguang.MLM_signals.csv").open("w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["SNP", "CHROM", "POS", "REF", "ALT", "MAF", "Effect", "SE", "jap_Ganguang.MLM"])
            writer.writerow(["snp1", "6", "8359017", "G", "A", "0.3", "3.0", "0.5", "1e-9"])
            writer.writerow(["snp2", "9", "16271846", "C", "T", "0.4", "2.2", "0.3", "1e-10"])

        subprocess.run([sys.executable, str(script), "jap_Ganguang"], cwd=root, check=True)

        rows = (root / "jap_Ganguang.MLM.peak").read_text(encoding="utf-8").splitlines()

    assert len(rows) == 2
    assert rows[0].split("\t")[:3] == ["snp1", "6", "8359017"]
    assert rows[1].split("\t")[:3] == ["snp2", "9", "16271846"]


def test_farmcpu_only_rap_intermediate_is_generated_from_peak_and_expression() -> None:
    script = REPO_ROOT / "V1" / "code" / "select_FarmCPUpeak_info_RAP.py"

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        basic = root / "basic_data"
        basic.mkdir()
        (root / "HNHZ.FarmCPU.peak").write_text(
            "Chromosome\tPosition\t-log10(Pvalue)\tParent_geno\t1_genotype\t2_genotype\tEffect\tSE\n"
            "1\t1200\t12.00\t*\tC\tT\t0.4\t0.04\n",
            encoding="utf-8",
        )
        (root / "HNHZ.winqtlcart.peak").write_text(
            "NAM1\t1\t0.001\t0.002\t4.5\n",
            encoding="utf-8",
        )
        lod_header = (
            "Geno_id\tChr\tGeno_star\tGeno_end\tNAM1\tNAM2\tNAM3\tNAM4\tNAM5\tNAM6\tNAM7\tNAM8\t"
            "NAM9\tNAM10\tNAM11\tNAM12\tNAM13\tNAM14\tNAM15\n"
        )
        (root / "HNHZ_NAM_LOD_geno.info").write_text(
            lod_header + "LOC_Os01g00010\t1\t1000\t1300\t4.5\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n",
            encoding="utf-8",
        )
        (root / "HNHZ_NAM_LOD_geno_rap.info").write_text(lod_header, encoding="utf-8")
        (basic / "combineRAPvsMSU_RNAseq.txt").write_text(
            "Gene_ID\tChr\tStart\tEnd\troot_seeding_stage\tshoot_seeding_stage\tleaf_sheath_seeding_stage\t"
            "pistils\tanthers\tLeaves_seeding_stage\tleaves_tillering_stage\tleaves_flowering_stage\t"
            "young_panicle\tpanicle_flowering_stage\tpanicle_filling_stage\ttiller_buds\tembryos\t"
            "developing_seeds\tseeds_after_ageing\tseed_Germinating_stage\tshoot_apical_meristem\tlamina_joint\n"
            "LOC_Os01g00010\t1\t1000\t1300\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\t13\t14\t15\t16\t17\t18\n",
            encoding="utf-8",
        )

        subprocess.run(
            [sys.executable, str(script), "HNHZ", "leaves_flowering_stage", "--candidate-window-bp", "500"],
            cwd=root,
            check=True,
        )

        rows = (root / "RAP_HNHZ.FarmCPUpeak_info").read_text(encoding="utf-8").splitlines()

    assert len(rows) == 2
    parts = rows[1].split("\t")
    assert parts[:7] == [
        "QTL1",
        "Chr1-1200",
        "12.00",
        "*",
        "NAM1:chr1:1000-2000:4.5;",
        "LOC_Os01g00010",
        "1:1000-1300",
    ]
    assert parts[7] == "4.5 0 0 0 0 0 0 0 0 0 0 0 0 0 0"
    assert parts[8:13] == ["*", "*", "*", "*", "*"]
    assert parts[17:] == [str(i) for i in range(1, 19)]


def test_farmcpu_only_rap_intermediate_fills_variant_scores() -> None:
    script = REPO_ROOT / "V1" / "code" / "select_FarmCPUpeak_info_RAP.py"

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        basic = root / "basic_data"
        basic.mkdir()
        (root / "HNHZ.FarmCPU.peak").write_text(
            "Chromosome\tPosition\t-log10(Pvalue)\tParent_geno\t1_genotype\t2_genotype\tEffect\tSE\n"
            "1\t1200\t12.00\t*\tC\tT\t0.4\t0.04\n",
            encoding="utf-8",
        )
        (root / "HNHZ.winqtlcart.peak").write_text(
            "NAM1\t1\t0.001\t0.002\t4.5\n",
            encoding="utf-8",
        )
        lod_header = (
            "Geno_id\tChr\tGeno_star\tGeno_end\tNAM1\tNAM2\tNAM3\tNAM4\tNAM5\tNAM6\tNAM7\tNAM8\t"
            "NAM9\tNAM10\tNAM11\tNAM12\tNAM13\tNAM14\tNAM15\n"
        )
        (root / "HNHZ_NAM_LOD_geno.info").write_text(
            lod_header + "LOC_Os01g00010\t1\t1000\t1300\t4.5\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n",
            encoding="utf-8",
        )
        (root / "HNHZ_NAM_LOD_geno_rap.info").write_text(lod_header, encoding="utf-8")
        (basic / "combineRAPvsMSU_RNAseq.txt").write_text(
            "Gene_ID\tChr\tStart\tEnd\troot_seeding_stage\tshoot_seeding_stage\tleaf_sheath_seeding_stage\t"
            "pistils\tanthers\tLeaves_seeding_stage\tleaves_tillering_stage\tleaves_flowering_stage\t"
            "young_panicle\tpanicle_flowering_stage\tpanicle_filling_stage\ttiller_buds\tembryos\t"
            "developing_seeds\tseeds_after_ageing\tseed_Germinating_stage\tshoot_apical_meristem\tlamina_joint\n"
            "LOC_Os01g00010\t1\t1000\t1300\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\t13\t14\t15\t16\t17\t18\n",
            encoding="utf-8",
        )
        sift_dir = root / "SIFTtolerantScore_MSUvsRAP.txt"
        sift_dir.mkdir()
        (sift_dir / "SIFTtolerantScore_MSUvsRAP.txt").write_text(
            "1\t1100\t1100\tC\tG\t0.01\tLOC_Os01g00010.1\tNONSYNONYMOUS\n",
            encoding="utf-8",
        )
        variant_dir = root / "21117007"
        variant_dir.mkdir()
        (variant_dir / "Chr1_combineSNPeffect_MSUvsRAP.txt").write_text(
            "#CHRO\tPos\tRef\tAlt\tEffect\tP1\tP2\n"
            "Chr1\t1100\tC\tG\t*\tC\tG\n",
            encoding="utf-8",
        )
        (variant_dir / "Chr1_combineINDELtolerantScore_MSUvsRAP.txt").write_text(
            "#CHRO\tPos_star\tPos_end\tRef\tAlt\tEffect\tP1\tP2\tScore\n"
            "Chr1\t1110\t1120\tA\t-\tLOC_Os01g00010\t0\t1\t17\n",
            encoding="utf-8",
        )
        (variant_dir / "NAM_Chr1combineRAPvsMSU_codingSV.txt").write_text(
            "GeneID\tChr\tPOS_star\tPOS_end\tvariation\tP1\tP2\tScore\n"
            "LOC_Os01g00010\t1\t1150\t1160\t10bp_de\t0|0\t1|1\t25\n",
            encoding="utf-8",
        )

        subprocess.run(
            [sys.executable, str(script), "HNHZ", "leaves_flowering_stage", "--candidate-window-bp", "500"],
            cwd=root,
            check=True,
        )

        rows = (root / "RAP_HNHZ.FarmCPUpeak_info").read_text(encoding="utf-8").splitlines()

    parts = rows[1].split("\t")
    assert parts[13:16] == ["20", "17", "25"]


def test_annotate_recommendations_adds_keywords_and_ranks_top_genes() -> None:
    script = REPO_ROOT / "V1" / "bin" / "annotate_recommendations.py"

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        basic = root / "basic_data"
        basic.mkdir()
        (basic / "all.locus_brief_info.7.0.with_keyword.tsv").write_text(
            "chr\tlocus\tmodel\tstart\tstop\tori\tis_TE\tis_expressed\tis_representative\tSymbol\tKeyword\tannotation\n"
            "Chr1\tLOC_Os01g00010\tLOC_Os01g00010.1\t100\t200\t+\tN\tY\tY\tHd1\tflowering;heading\tHeading date regulator\n"
            "Chr1\tLOC_Os01g00020\tLOC_Os01g00020.1\t300\t400\t+\tN\tY\tY\t\t\tExpressed protein\n",
            encoding="utf-8",
        )
        (root / "win.tsv").write_text(
            "QTL\tFarmCPU_pos\tFarmCPU_-logP\tMLM_INFO(Chr-pos-(-log10))\tGeno\tGeno_pos\t"
            "SIFT_Score\tINDEL_Score\tSV_Score\tScore_peak\tScore_TE\tScore_annotate\tScore_expression\tScore_match\tScore\n"
            "QTL1\tChr1-150\t12.0\t*\tLOC_Os01g00010\t1:100-200\t25\t10\t5\t10\t0\t20\t15\t5\t60\n"
            "QTL1\tChr1-350\t10.0\t*\tLOC_Os01g00020\t1:300-400\t1\t*\t*\t10\t0\t0\t5\t*\t20\n",
            encoding="utf-8",
        )
        (root / "window.tsv").write_text(
            "QTL\tFarmCPU_pos\tFarmCPU_-logP\tMLM_INFO(Chr-pos-(-log10))\tGeno\tGeno_pos\t"
            "SIFT_Score\tINDEL_Score\tSV_Score\tScore_peak\tScore_TE\tScore_annotate\tScore_expression\tScore_match\tScore\n"
            "QTL1\tChr1-150\t12.0\t*\tLOC_Os01g00010\t1:100-200\t20\t10\t5\t8\t0\t20\t10\t*\t50\n",
            encoding="utf-8",
        )

        subprocess.run(
            [
                sys.executable,
                str(script),
                "--inputs",
                "win.tsv",
                "window.tsv",
                "--keyword-table",
                str(basic / "all.locus_brief_info.7.0.with_keyword.tsv"),
                "--output-dir",
                str(root),
                "--output-prefix",
                "HNHZ",
            ],
            cwd=root,
            check=True,
        )

        annotated = (root / "win.annotated.tsv").read_text(encoding="utf-8").splitlines()
        top = (root / "HNHZ_top10_candidate_genes.tsv").read_text(encoding="utf-8").splitlines()

    annotated_header = annotated[0].split("\t")
    annotated_row = annotated[1].split("\t")
    assert annotated_row[annotated_header.index("Symbol")] == "Hd1"
    assert annotated_row[annotated_header.index("Keyword")] == "flowering;heading"
    assert annotated_row[annotated_header.index("Locus_annotation")] == "Heading date regulator"

    top_header = top[0].split("\t")
    top_row = top[1].split("\t")
    assert top_row[top_header.index("Gene")] == "LOC_Os01g00010"
    assert top_row[top_header.index("Modes")] == "win;window"
    assert top_row[top_header.index("Symbol")] == "Hd1"


def test_step2_tgw_keyword_uses_grain_weight_terms_and_developing_seed_expression() -> None:
    script = REPO_ROOT / "V1" / "code" / "RAP_Stpe2_select_FarmCPUpeak_info.py"

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        basic = root / "basic_data"
        basic.mkdir()
        (basic / "key_word.txt").write_text(
            "Yield:grain,seed,weight;\n",
            encoding="utf-8",
        )

        fields = [
            "QTL1",
            "Chr1-1000",
            "10.0",
            "*",
            "*",
            "LOC_Os01g00010",
            "1:100-200",
            "*",
            "*",
            "*",
            "*",
            "*",
            "grain weight regulator",
            "*",
            "*",
            "*",
            "*",
        ]
        expression = ["0.001"] * 18
        expression[13] = "1.0"  # developing_seeds
        (root / "RAP_all_TGW_quality.FarmCPUpeak_info").write_text(
            "header\n" + "\t".join(fields + expression) + "\n",
            encoding="utf-8",
        )

        subprocess.run([sys.executable, str(script), "all_TGW_quality", "TGW"], cwd=root, check=True)

        rows = (root / "RAP_Step2_all_TGW_quality.FarmCPUpeak_info").read_text(encoding="utf-8").splitlines()

    header = rows[0].split("\t")
    row = rows[1].split("\t")
    assert float(row[header.index("Score_annotate")]) > 0
    assert float(row[header.index("Score_expression")]) > 0


def test_step2_heading_alias_uses_heading_date_terms_and_flowering_expression() -> None:
    script = REPO_ROOT / "V1" / "code" / "RAP_Stpe2_select_FarmCPUpeak_info.py"

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        basic = root / "basic_data"
        basic.mkdir()
        (basic / "key_word.txt").write_text(
            "Heading_date:flowering,heading date,photoperiod;\n",
            encoding="utf-8",
        )

        fields = [
            "QTL1",
            "Chr6-8359017",
            "36.0",
            "*",
            "*",
            "LOC_Os06g00010",
            "6:100-200",
            "*",
            "*",
            "*",
            "*",
            "*",
            "flowering heading date regulator",
            "*",
            "*",
            "*",
            "*",
        ]
        expression = ["0.001"] * 18
        expression[7] = "1.0"  # leaves_flowering_stage
        (root / "RAP_jap_Ganguang.FarmCPUpeak_info").write_text(
            "header\n" + "\t".join(fields + expression) + "\n",
            encoding="utf-8",
        )

        subprocess.run([sys.executable, str(script), "jap_Ganguang", "heading"], cwd=root, check=True)

        rows = (root / "RAP_Step2_jap_Ganguang.FarmCPUpeak_info").read_text(encoding="utf-8").splitlines()

    header = rows[0].split("\t")
    row = rows[1].split("\t")
    assert float(row[header.index("Score_annotate")]) > 0
    assert float(row[header.index("Score_expression")]) > 0


if __name__ == "__main__":
    test_rmvp_signals_are_converted_to_farmcpu_peak()
    test_data_tgw_farmcpu_csv_is_filtered_to_default_significant_sites()
    test_mlm_peak_uses_named_pvalue_column_when_maf_is_present()
    test_mlm_signals_from_data_flowering_are_supported()
    test_farmcpu_only_rap_intermediate_is_generated_from_peak_and_expression()
    test_farmcpu_only_rap_intermediate_fills_variant_scores()
    test_annotate_recommendations_adds_keywords_and_ranks_top_genes()
    test_step2_tgw_keyword_uses_grain_weight_terms_and_developing_seed_expression()
    test_step2_heading_alias_uses_heading_date_terms_and_flowering_expression()
