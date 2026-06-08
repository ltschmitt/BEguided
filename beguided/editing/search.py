"""PAM search operations."""


def find_pams_regex(seq, start, end, pam_regex, pam_length):
    """Find all PAM start positions in search range using regex with lookahead.

    Args:
        seq: Genome sequence string.
        start: 0-based start of search range.
        end: 0-based end of search range.
        pam_regex: Compiled regex pattern.
        pam_length: Length of PAM sequence.

    Returns:
        list of 0-based PAM start positions.
    """
    pams = []
    search_seq = seq[start:end]
    seq_len = len(seq)

    for match in pam_regex.finditer(search_seq):
        # With lookahead pattern, match.start(1) is actual PAM position
        pam_start = start + match.start(1)
        if pam_start + pam_length <= seq_len:
            pams.append(pam_start)

    return pams
