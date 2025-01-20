import os
import re


def count_and_assign_go_annotations(go2protein_count, go_annotations, num_of_proteins):
    for go_id in go_annotations:
        count = go2protein_count.setdefault(go_id, 0)
        count += 1 * num_of_proteins
        go2protein_count[go_id] = count


def parse_ips_file(ips_file):

    go2protein_count = {}
    num_of_proteins_with_go = 0
    total_num_of_proteins = 0
    if os.path.exists(ips_file):
        handle = open(ips_file, "r")
        go_pattern = re.compile("GO:\\d+")
        line_counter = 0
        previous_protein_acc = None
        go_annotations_single_protein = set()
        for line in handle:
            line_counter += 1
            line = line.strip()
            chunks = line.split("\t")
            # Get protein accession
            current_protein_acc = chunks[0]
            num_of_proteins = len(current_protein_acc.split("|"))
            # If new protein accession extracted, store GO annotation counts in result dictionary
            if not current_protein_acc == previous_protein_acc:
                total_num_of_proteins += 1
                if len(go_annotations_single_protein) > 0:
                    num_of_proteins_with_go += 1

                previous_protein_acc = current_protein_acc
                count_and_assign_go_annotations(
                    go2protein_count, go_annotations_single_protein, num_of_proteins
                )
                # reset go id set because we hit a new protein accession
                go_annotations_single_protein = set()
            # Parse out GO annotations
            # GO annotations are associated to InterPro entries (InterPro entries start with 'IPR')
            # Than use the regex to extract the GO Ids (e.g. GO:0009842)
            if len(chunks) >= 13 and chunks[11].startswith("IPR"):
                for go_annotation in go_pattern.findall(line):
                    go_annotations_single_protein.add(go_annotation)

        # Do final counting for the last protein
        count_and_assign_go_annotations(
            go2protein_count, go_annotations_single_protein, num_of_proteins
        )
        total_num_of_proteins += 1

        handle.close()

    return go2protein_count
