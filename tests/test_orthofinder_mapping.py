import pandas as pd

from python.orthology.build_id_symbol_map_from_gtf import (
    build_mapping,
)
from python.orthology.map_orthofinder_ids_to_symbols import (
    annotate_cell,
    annotate_orthogroups,
)


def test_build_mapping_from_gtf(tmp_path):
    gtf = tmp_path / "test.gtf"
    gtf.write_text(
        'chr1\tsrc\tgene\t1\t10\t.\t+\t.\tgene_id "G1"; gene_name "GENE1";\n'
        'chr1\tsrc\tgene\t20\t30\t.\t+\t.\tgene_id "G2"; gene_name "GENE2";\n'
    )

    mapping = build_mapping(gtf)

    assert mapping.to_dict("records") == [
        {"id": "G1", "symbol": "GENE1"},
        {"id": "G2", "symbol": "GENE2"},
    ]


def test_annotate_cell_preserves_unmapped_ids():
    annotated, mapped, unmapped = annotate_cell(
        "G1, G2",
        {"G1": "GENE1"},
    )

    assert annotated == "GENE1|G1, UNMAPPED|G2"
    assert mapped == 1
    assert unmapped == 1


def test_annotate_orthogroups_uses_species_mapping_files(tmp_path):
    mapping_dir = tmp_path / "mappings"
    mapping_dir.mkdir()

    pd.DataFrame({
        "id": ["H1"],
        "symbol": ["GENE1"],
    }).to_csv(mapping_dir / "human.tsv", sep="\t", index=False)

    pd.DataFrame({
        "id": ["M1"],
        "symbol": ["Gene1"],
    }).to_csv(mapping_dir / "mouse.tsv", sep="\t", index=False)

    orthogroups = pd.DataFrame({
        "Orthogroup": ["OG1"],
        "human": ["H1"],
        "mouse": ["M1"],
    })

    annotated, summary = annotate_orthogroups(
        orthogroups,
        mapping_dir,
    )

    assert annotated.loc[0, "human"] == "GENE1|H1"
    assert annotated.loc[0, "mouse"] == "Gene1|M1"
    assert summary["unmapped_ids"].sum() == 0
