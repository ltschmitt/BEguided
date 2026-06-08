"""CRISPRater on-target score calculator.

Code adapted from CRISPRater (CCTOP standalone):
https://bitbucket.org/juanlmateo/cctop_standalone/src/master/cctop/scoreCRISPRater.py
"""

import re


MODEL_WEIGHT = [
    0.14177385,
    0.06966514,
    0.04216254,
    0.03303432,
    0.02355430,
    -0.04746424,
    -0.04878001,
    -0.06981921,
    -0.07087756,
    -0.08160700,
]
MODEL_OFFSET = 0.6505037

_PATTERN_CG = re.compile("[CG]")


def _get_gc_freq(seq):
    """Return GC content of *seq* as a float."""
    cg = len(_PATTERN_CG.findall(seq))
    return float(cg) / len(seq)


def _calc_features(seq):
    """Extract the 10 binary / continuous features from a 20-nt sequence."""
    if len(seq) != 20:
        raise ValueError(
            f"CRISPRater expects a 20-nt sequence, got {len(seq)}"
        )
    feat = [0] * 10
    feat[0] = _get_gc_freq(seq[3:13])
    if seq[19] == "G":
        feat[1] = 1
    if seq[2] in ("T", "A"):
        feat[2] = 1
    if seq[11] in ("G", "A"):
        feat[3] = 1
    if seq[5] == "G":
        feat[4] = 1
    if seq[3] in ("T", "A"):
        feat[5] = 1
    if seq[17] in ("G", "A"):
        feat[6] = 1
    if seq[4] in ("C", "A"):
        feat[7] = 1
    if seq[13] == "G":
        feat[8] = 1
    if seq[14] == "A":
        feat[9] = 1
    return feat


def get_score(seq):
    """Calculate CRISPRater on-target score for a 20-nt guide sequence.

    Args:
        seq: DNA sequence, uppercase, exactly 20 nucleotides.

    Returns:
        float: Predicted on-target activity score.

    Raises:
        ValueError: If *seq* is not exactly 20 nt.
    """
    features = _calc_features(seq)
    score = sum(f * w for f, w in zip(features, MODEL_WEIGHT)) + MODEL_OFFSET
    return score
