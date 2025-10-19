import edlib

def find_substring_with_mismatches(dna_string, substring, max_mismatches):
    # Perform the alignment
    result = edlib.align(substring, dna_string, mode="HW", task="locations", k=max_mismatches)

    if result['locations']:
        start, end = result['locations'][0]
        return start, end + 1  # edlib returns inclusive end index

    return None, None

# Example usage
dna_string = "ACTGATGATGTGGACCTGGAACAGGTGCGGCAGCTGGTGCCTCG"
substring = "ATGATGTGGACCTGGAACAG"
max_mismatches = 2

upstream_exon = 'ATGATGTGGACCTGGAACAG'  # last 20 bp
downstream_exon = 'GTGCGGCAGCTGGTGCCTCG' # first 20 bp

start, end = find_substring_with_mismatches(dna_string, substring, max_mismatches)

if start is not None and end is not None:
    print(f"Substring found from index {start} to {end-1}")
else:
    print("Substring not found within the allowed mismatches")