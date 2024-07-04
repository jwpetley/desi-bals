rule fractions:
    input:
        "src/data/desi_bals_radio.fits"
        "src/data/desi_qso_radio.fits"
    output:
        "src/tex/figures/detection_fraction_z.pdf"
    script:
        "scr/scripts/detection_fraction.py"


rule matching:
    output:
        directory("src/data/matched")
    cache:
        True
    script:
        "scr/scripts/matching.py"

rule matching_figures:
    input:
        "src/data/desi_bals_radio.fits"
        "src/data/desi_qso_radio.fits"
    output:
        "src/tex/figures/matched_detection_fraction.pdf"
        "src/tex/figures/six_micron_z.pdf"
    script:
        "scr/scripts/matching.py"

rule correlation:
    input:
        "src/data/matched"
    output:
        "scr/figures/flux_bi.pdf"
    script:
        "scr/scripts/correlations.py"