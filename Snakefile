rule matching:
    output:
        directory("src/data/matched")
    cache:
        True
    script:
        "scr/scripts/matching.py"

rule correlation:
    input:
        "src/data/matched"
    output:
        "scr/figures/flux_bi.pdf"
    script:
        "scr/scripts/correlations.py"