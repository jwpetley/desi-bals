rule fractions:
    input:
        "src/data/desi_bals_radio.fits"
        "src/data/desi_qso_radio.fits"
    output:
        "src/tex/figures/detection_fraction_z.pdf"
    script:
        "src/scripts/detection_fraction.py"


rule matching:
    output:
        directory("src/data/matched")
    cache:
        True
    script:
        "src/scripts/matching.py"

rule matching_figures:
    input:
        "src/data/matched"
        "src/data/desi_bals_radio.fits"
        "src/data/desi_qso_radio.fits"
    output:
        "src/tex/figures/matched_detection_fraction.pdf"
        "src/tex/figures/six_micron_z.pdf"
    script:
        "src/scripts/matching_figures.py"

rule correlation:
    input:
        "src/data/matched"
    output:
        "src/figures/flux_bi.pdf"
    script:
        "src/scripts/correlations.py"