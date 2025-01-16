import argparse
import os
import re

from go_utils import parse_ips_file

description = "Go slim pipeline."
parser = argparse.ArgumentParser(description=description)
parser.add_argument("-go", "--go_obo", help="Gene Ontology basic file.", required=True)
parser.add_argument(
    "-gb", "--go_banding", help="Subset GO banding file.", required=True
)
parser.add_argument(
    "-i", "--ips_input", help="InterProScan result file.", required=True
)
parser.add_argument("-o", "--output", help="GO summary output file.", required=True)
args = parser.parse_args()

GO_OBO = args.go_obo
GO_BANDING = args.go_banding
IPS_INPUT = args.ips_input
OUTPUT = args.output


def parse_mapped_gaf_file(gaf_file):
    """
    parse_mapped_gaf_file(gaf_file) -> dictionary
    Example of GAF mapped output:
        !gaf-version: 2.0
        ! This GAF has been mapped to a subset:
        ! Subset: user supplied list, size = 38
        ! Number of annotation in input set: 1326
        ! Number of annotations rewritten: 120
        EMG	GO:0005839	GO		GO:0005839	PMID:12069591	IEA		C			protein	taxon:1310605	20160528	InterPro
        EMG	GO:0000160	GO		GO:0005575	PMID:12069591	IEA		C			protein	taxon:1310605	20160528	InterPro
    Parsing the above GAF file will create the following dictionary:
    result = {'GO:0005839':'GO:0005839', 'GO:0000160':'GO:0005575'}
    :param gaf_file:
    :return:
    """
    result = {}
    if os.path.exists(gaf_file):
        handle = open(gaf_file, "r")
        for line in handle:
            if not line.startswith("!"):
                line = line.strip()
                splitted_line = line.split("\t")
                go_id = splitted_line[1]
                mapped_go_id = splitted_line[4]
                result.setdefault(go_id, set()).add(mapped_go_id)
    return result


def parse_iprscan_output_goslim_counts(iprscanOutput, map2slim_mapped_go_ids_dict):
    # result -> GO accessions mapped to number of occurrences
    # Example {'GO:0009842':267, 'GO:0009841':566}
    result = {}
    if os.path.exists(iprscanOutput):
        handle = open(iprscanOutput, "r")
        # Example GO Id -> GO:0009842
        goPattern = re.compile("GO:\\d+")
        line_counter = 0
        previous_protein_acc = None
        go_annotations_single_protein = set()
        # Set default value for number of proteins to 1
        num_of_proteins = 1
        for line in handle:
            line_counter += 1
            line = line.strip()
            chunks = line.split("\t")
            # Get protein accession
            current_protein_acc = chunks[0]
            num_of_proteins = len(current_protein_acc.split("|"))
            # If new protein accession extracted, store GO annotation counts in result dictionary
            if not current_protein_acc == previous_protein_acc:
                previous_protein_acc = current_protein_acc
                count_slims(
                    go_annotations_single_protein,
                    map2slim_mapped_go_ids_dict,
                    num_of_proteins,
                    result,
                )
                # reset go id set because we hit a new protein accession
                go_annotations_single_protein = set()
            # Parse out GO annotations
            # GO annotations are associated to InterPro entries (InterPro entries start with 'IPR')
            # Than use the regex to extract the GO Ids (e.g. GO:0009842)
            if len(chunks) >= 13 and chunks[11].startswith("IPR"):
                for go_annotation in goPattern.findall(line):
                    go_annotations_single_protein.add(go_annotation)

        # Do final counting for the last protein
        count_slims(
            go_annotations_single_protein,
            map2slim_mapped_go_ids_dict,
            num_of_proteins,
            result,
        )
        handle.close()
    return result


def count_slims(
    go_annotations_single_protein, map2slim_mapped_go_ids_dict, num_of_proteins, result
):
    # count goslims
    slim_go_ids_set = set()
    # Get the set of slim terms
    for go_annotation in go_annotations_single_protein:
        mapped_go_ids = map2slim_mapped_go_ids_dict.get(go_annotation)
        if mapped_go_ids:
            slim_go_ids_set.update(mapped_go_ids)
    # Iterate over the set of slim terms and update the counts
    for slim_go_id in slim_go_ids_set:
        count = result.setdefault(slim_go_id, 0)
        count += 1 * num_of_proteins
        result[slim_go_id] = count


def get_go_slim_summary(go_slim_banding_file, go_slims_2_protein_count):
    summary = []

    file_handler = open(go_slim_banding_file, "r")

    for line in file_handler:
        if line.startswith("GO"):
            line = line.strip()
            line_chunks = line.split("\t")
            go_id = line_chunks[0]
            term = line_chunks[1]
            category = line_chunks[2]
            # Default value for the count
            count = 0
            if go_id in go_slims_2_protein_count:
                count = go_slims_2_protein_count.get(go_id)
            summary.append((go_id, term, category, count))
    return summary


def write_go_summary_to_file(goSummary, outputFile):
    handle = open(outputFile, "w")
    for go, term, category, count in goSummary:
        handle.write('","'.join(['"' + go, term, category, str(count) + '"']) + "\n")
    handle.close()


def get_gene_ontology(obo_file):
    """
    Parses OBO formatted file.
    :param obo_file:
    :return:
    """
    result = []
    handle = open(obo_file, "r")
    id, term, category = "", "", ""
    for line in handle:
        line = line.strip()
        splitLine = line.split(": ")
        if line.startswith("id:"):
            id = splitLine[1].strip()
        elif line.startswith("name:"):
            term = splitLine[1].strip()
        elif line.startswith("namespace"):
            category = splitLine[1].strip()
        else:
            if id.startswith("GO:") and id and term and category:
                item = (id, term, category)
                result.append(item)
                id, term, category = "", "", ""
    handle.close()
    return result


def go_sort_key(item):
    return (item[2], -item[3])


def get_full_go_summary(core_gene_ontology, go2protein_count_dict, topLevelGoIds):
    summary = []

    for goId, term, category in core_gene_ontology:

        if (goId in go2protein_count_dict) and (
            goId not in topLevelGoIds
        ):  # make sure that top level terms are not included (they tell you nothing!)
            count = go2protein_count_dict.get(goId)
            summary.append((goId, term, category, count))
    summary.sort(key=go_sort_key)
    return summary


def main():

    print("Parsing the InterProScan result output file: " + IPS_INPUT)
    go2protein_count_dict = parse_ips_file(IPS_INPUT)
    print("Finished parsing.")

    # Generate GO summary
    print("Loading full Gene ontology: " + GO_OBO)
    core_gene_ontology_list = get_gene_ontology(GO_OBO)
    print("Finished loading.")

    print("Generating full GO summary...")
    topLevelGoIds = ["GO:0008150", "GO:0003674", "GO:0005575"]
    full_go_summary = get_full_go_summary(
        core_gene_ontology_list, go2protein_count_dict, topLevelGoIds
    )
    # delete core gene ontology list
    del core_gene_ontology_list
    print("Finished generation.")

    print("Writing full GO summary to the following file: " + OUTPUT)
    write_go_summary_to_file(full_go_summary, OUTPUT)
    # delete full GO summary variable
    del full_go_summary
    print("Finished writing.")

    go2mapped_go = parse_mapped_gaf_file("ERR435123_goslim_annotations.gaf")
    print("Finished parsing.")

    print("Getting GO slim counts by parsing I5 TSV again")
    go_slims_2_protein_count = parse_iprscan_output_goslim_counts(
        IPS_INPUT, go2mapped_go
    )

    go_slim_summary = get_go_slim_summary(GO_BANDING, go_slims_2_protein_count)
    go_slim_output_file = OUTPUT + "_slim"
    print("Writing GO slim summary to the following file: " + go_slim_output_file)
    write_go_summary_to_file(go_slim_summary, go_slim_output_file)
    # delete full GO summary variable
    del go_slim_summary
    print("Finished writing.")


if __name__ == "__main__":
    main()
