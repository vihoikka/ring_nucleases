import pandas as pd
import re
import openpyxl

# CARF family motifs dicitonary based on Makarova 2020 Figure 4.
motifs = {
    "CARF7_Motif-I": r"GTS",
    "CARF_m4_Motif-I": r"G[LV][ST]",
    "CARF5_Motif-I": r"G[LAM][ST]",
    "CARF_m13_Motif-I": r"GTSP",
    "CARF1a_Motif-I": r"G[TG][ST]",
    "CARF1b_Motif-I": r"G[NTD][TS]",
    "CARF9_Motif-I": r"G[ND]F",
    "CARF7_Motif-IA": r"SAE",
    "CARF_m4_Motif-IA": r"[PA]ELH",
    "CARF7_Motif-II": r"N[ALP][TY]GG[FY]K",
    "CARF_m4_Motif-II": r"[SG][ILV]T[ST]G[FAW][NKR]",
    "CARF5_Motif-II": r"S[IL]AGGRK",
    "CARF_m13_Motif-II": r"D[IVF]TGGRK",
    "CARF1a_Motif-II": r"D.T[GS]GTK",
    "CARF1b_Motif-II": r"[NS]...GT",
    "CARF9_Motif-II": r"D[LST]THGLN",
}

# Grouping motifs into CARF families
carf_families = {
    "CARF7": ["CARF7_Motif-I", "CARF7_Motif-IA", "CARF7_Motif-II"],
    "CARF_m4": ["CARF_m4_Motif-I", "CARF_m4_Motif-IA", "CARF_m4_Motif-II"],
    "CARF5": ["CARF5_Motif-I", "CARF5_Motif-II"],
    "CARF_m13": ["CARF_m13_Motif-I", "CARF_m13_Motif-II"],
    "CARF1a": ["CARF1a_Motif-I", "CARF1a_Motif-II"],
    "CARF1b": ["CARF1b_Motif-I", "CARF1b_Motif-II"],
    "CARF9": ["CARF9_Motif-I", "CARF9_Motif-II"],
}


# Function to check for motifs in a sequence
def check_motifs(sequence, motifs):
    results = {}
    for motif_name, motif_pattern in motifs.items():
        # Use regex to check if the motif is present in the sequence
        if re.search(motif_pattern, sequence):
            results[motif_name] = True
        else:
            results[motif_name] = False
    return results


# Function to determine the CARF family based on motif support
def determine_carf_family(row, carf_families):
    supported_families = []
    for family, family_motifs in carf_families.items():
        # Check if any motif in the family is supported
        if any(row[motif] for motif in family_motifs if motif in row):
            supported_families.append(family)

    # Determine the consensus
    if len(supported_families) == 1:
        return supported_families[0]  # Single family match
    elif len(supported_families) > 1:
        # investigate if one family has more support than any other
        family_support = {family: 0 for family in supported_families}
        for family in supported_families:
            for motif in carf_families[family]:
                if row[motif]:
                    family_support[family] += 1
        max_support = max(family_support.values())
        if list(family_support.values()).count(max_support) == 1:
            return max(family_support, key=family_support.get)
        else:
            return "No consensus"
    else:
        return "No matches"  # No motifs supported


# Load the Excel file into a pandas DataFrame
def load_excel(file_path):
    return pd.read_excel(file_path)


# Main function to process the sequences and check motifs
def process_sequences(file_path, output_path):
    # Load the Excel file
    df = load_excel(file_path)

    # Ensure the column "RN_sequence" exists
    if "RN_sequence" not in df.columns:
        raise ValueError("The input file must contain a column named 'RN_sequence'")

    # Create columns for each motif
    for motif_name in motifs.keys():
        df[motif_name] = df["RN_sequence"].apply(
            lambda seq: check_motifs(seq, motifs)[motif_name]
        )

    print("Motifs checked!")

    # Add a column for the CARF family consensus
    df["CARF_Family"] = df.apply(
        determine_carf_family, axis=1, carf_families=carf_families
    )

    # Save the results to a new Excel file
    df.to_excel(output_path, index=False)
    print(f"Results saved to {output_path}")


if __name__ == "__main__":
    # Input Excel file path
    input_file = "supplementary_rn_file.xlsx"
    output_file = "supplementary_rn_file_with_motifs_and_consensus.xlsx"

    # Process the sequences and check for motifs
    process_sequences(input_file, output_file)
