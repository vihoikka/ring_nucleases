from BCBio import GFF
from dna_features_viewer import GraphicFeature, GraphicRecord, BiopythonTranslator
import matplotlib.pyplot as plt
import argparse
import pandas as pd
import os
import gffutils

"""
A modified version of locus_visualiser.py that draws genomic maps of RNs and their neighbouring genes outside of CRISPR loci.
Only uses the genomic gff file to draw the maps
"""

# create args inputs for cas_operons_file, locus, output_folder, host_genomes_folder, mode, hmm_rows, hmm_targets, catyper_out, cctyper_protein_table
parser = argparse.ArgumentParser(description="Visualises expanded CRISPR loci")
parser.add_argument("-rnw", "--rn_wildcard", help="input file", required=True)
parser.add_argument("-rnf", "--rn_file", help="rn file tsv", required=True)
parser.add_argument("-o", "--output_folder", help="output folder", required=True)
parser.add_argument(
    "-hg", "--host_genomes_folder", help="host genomes folder", required=True
)

args = parser.parse_args()
print("Starting RN viz script")

# generate bash command for the above
# python locus_visualiser.py -i cas_operons.tsv -l locus -o output_folder -hg host_genomes_folder

# assign args to similarly named variables
rn_wildcard = args.rn_wildcard
rn_file = args.rn_file
output_folder = args.output_folder
host_genomes_folder = args.host_genomes_folder

rn_df = pd.read_csv(rn_file, sep="\t")

# when looking for effectors in the locus, the effector search range is the number of bases up or downstream of the RN
search_range = 20000

# get sample name from first row of rn["host"]
sample = rn_df["host"][0]

gff_file = os.path.join(host_genomes_folder, sample, sample + "_features.gff")

# get contig ID
contig = rn_df["diamond_contig"][0]
print("Contig: " + contig)

# define start and end coordinates for plotting
viz_start = int(rn_df["genomic_start"][0]) - search_range
viz_end = int(rn_df["genomic_end"][0]) + search_range
features_list = []

print("Viz start coordinate: " + str(viz_start))
print("Viz end coordinate: " + str(viz_end))
print("RN start coordinate: " + str(rn_df["genomic_start"][0]))

rn_start = int(rn_df["genomic_start"][0])


def create_gff_iterator(gff_file, contig):
    """
    Creates a GFF iterator for a given GFF file and contig.
    Note that the if the in_handle is closed, then the script will crash later when the iterator is used. This is why we leave it open
    """
    in_handle = open(gff_file)
    limit_info = dict(
        gff_id=[contig], gff_type=["CDS", "repeat_region"]
    )  # Add or remove types as needed
    gff_iterator = GFF.parse(in_handle, limit_info=limit_info)
    return gff_iterator


def add_features_from_gff(
    gff_iterator, features_list, cas_operon_start, cas_operon_end, source, rn_start
):
    """
    Adds GraphicFeature objects to a list from a given GFF iterator
    """
    color_dict = {"original_gff": "#ffcccc", "rn": "#ffe747"}
    new_features_list = []
    for rec in gff_iterator:
        for feature in rec.features:
            if (
                int(feature.location.start) >= cas_operon_start
                and int(feature.location.end) <= cas_operon_end
            ):
                print(
                    "Feature with start: "
                    + str(feature.location.start)
                    + " and end: "
                    + str(feature.location.end)
                )
                label = feature.qualifiers.get("product", [None])[0]
                # if label is not NoneType and contains the string "unknown" or "hypothetical", then make label None
                if label is not None and (
                    "unknown" in label or "hypothetical" in label
                ):
                    label = None
                true_rn_start = rn_start - 1
                if int(feature.location.start) == true_rn_start:
                    print("Feature is the original RN")
                    color = color_dict["rn"]
                else:
                    color = color_dict[source]
                strand = 1 if feature.strand == 1 else -1
                start_coord_orf = int(feature.location.start) - cas_operon_start
                end_coord_orf = int(feature.location.end) - cas_operon_start
                print("Start coord ORF: " + str(start_coord_orf))
                print("End coord ORF: " + str(end_coord_orf))
                print("Strand: " + str(strand))
                print("Name: " + str(label))
                print("----------------")
                graphic_feature = GraphicFeature(
                    start=start_coord_orf,
                    end=end_coord_orf,
                    strand=strand,
                    color=color,
                    label=label,
                )
                new_features_list.append(graphic_feature)
    # combine new and old feature lists
    features_list = features_list + new_features_list
    return features_list


gff_iterator = create_gff_iterator(gff_file, contig)

# Prepare a list to hold the Graphic Features
features_list = []

# Add features from GFF
features_list = add_features_from_gff(
    gff_iterator,
    features_list,
    viz_start,
    viz_end,
    source="original_gff",
    rn_start=rn_start,
)

print("Number of features: " + str(len(features_list)))
# print all feature products

# Generate plot using DNA Features Viewer
record = GraphicRecord(sequence_length=viz_end - viz_start, features=features_list)
ax, _ = record.plot(figure_width=10)

# save
plt.savefig(output_folder + "/" + rn_wildcard + "_rn_viz.png")

# make a copy with info on whether this is on a crispr locus embedded in the filename
print("CRISPR status: " + str(rn_df["within_crispr"][0]))
if rn_df["within_crispr"][0]:
    plt.savefig(output_folder + "/" + rn_wildcard + "_rn_viz_crisprTrue.png")
else:
    plt.savefig(output_folder + "/" + rn_wildcard + "_rn_viz_crisprFalse.png")
