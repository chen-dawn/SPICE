import numpy as np

IN_MAP = np.asarray(
    [[0, 0, 0, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
)


def get_rev_comp(dna_sequence):
    """Generate the reverse complement of a DNA sequence.

    Args:
        dna_sequence (str): A string representing the DNA sequence (only A, T, C, and G).

    Returns:
        str: The reverse complement of the DNA sequence.
    """
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}

    if any(nucleotide not in complement for nucleotide in dna_sequence):
        raise ValueError(
            "Invalid DNA sequence. Sequence must contain only A, T, C, and G."
        )

    reverse_comp = "".join(
        complement[nucleotide] for nucleotide in reversed(dna_sequence)
    )

    return reverse_comp


def one_hot_encode(Xd):
    """
    One-hot encode the sequence.
    """
    return IN_MAP[Xd.astype("int8")]


def create_datapoint(seq, celltype, skipped_count, included_count):
    """
    Create a datapoint from a sequence, celltype, skipped_count and included_count.
    """
    seq = seq.decode("utf-8")
    seq = seq.upper().replace("A", "1").replace("C", "2")
    seq = seq.replace("G", "3").replace("T", "4").replace("N", "0")
    X0 = np.array(list(map(int, list(seq))))
    one_hot_encode_seq = one_hot_encode(X0)
    Yd = np.array([skipped_count, included_count]).transpose()
    return np.array([one_hot_encode_seq]), Yd.reshape(1, 2)  # X and Y.
