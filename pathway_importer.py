from gocamgen.gocamgen import GoCamModel
from ontobio.vocabulary.relations import OboRO
from rdflib.term import URIRef
from rdflib.namespace import Namespace, OWL
from prefixcommons.curie_util import expand_uri
from pathway_connections import PathwayConnectionSet
import argparse

ro = OboRO()
ENABLED_BY = URIRef(expand_uri(ro.enabled_by))

parser = argparse.ArgumentParser()
parser.add_argument('-f', "--filename", type=str, required=True,
                    help="Input filename of Signor pathway data")

def model_contains_statement(model, subject_uri, rel, object_id):
    for uri in model.uri_list_for_individual(object_id):
        if object_id == "GO:0004672":
            print(uri)
        axiom = find_statement_axiom(model, (subject_uri,rel,uri))
        # if model.writer.writer.graph.__contains__((subject_uri,rel,uri)):
        if axiom is not None:
            return True
    return False

def find_statement_axiom(model, triple):
    found_triples = model.writer.writer.graph.triples(triple)
    for found_one in found_triples:
        return found_one

def test_label_finding(model):
    # target = "UniProtKB:P84022"  # SMAD3
    target = "UniProtKB:P01106"  # MYC
    axiom_counter = 1
    for axiom in model.axioms_for_source(target, OWL.annotatedTarget):
        print("Axiom " + str(axiom_counter))
        for t in model.writer.writer.graph.triples((axiom, None, None)):    # t=(axiom,property,ind)
            labels = model.individual_label_for_uri(t[2])
            if len(labels) > 0:
                print(t[2])
                print(t[1])
                print(labels)
        axiom_counter += 1


def main():
    args = parser.parse_args()

    model = GoCamModel("superfamily_test.ttl")
    # p_connections = PathwayConnectionSet("SIGNOR-G2-M_trans_02_03_18.tsv")
    p_connections = PathwayConnectionSet(args.filename)
    linenum = 1
    # complex_csv_filename = "SIGNOR_complexes.csv"
    # complexes = SignorComplexFactory(complex_csv_filename).complexes

    total_pcs = len(p_connections.connections)
    print(total_pcs)
    skipped_count = 0

    # fill in regulated activities
    for pc in p_connections.connections:
        # if pc.id_a.startswith("SIGNOR") or pc.id_b.startswith("SIGNOR"):
        #     # for now to see how model first looks - skip complexes
        #     continue
        regulated_activity_pc = p_connections.find_other_regulated_activity(pc.id_b)  # find_by_id_a(p_connections.connections, pc.id_b) - regulated_activity_pc.mechanism["term"]
        if regulated_activity_pc is not None:
            regulated_activity_term = regulated_activity_pc.mechanism["term"]
            # regulated_activity_term_uri = regulated_activity_pc.individuals[regulated_activity_pc.mechanism["term"]]
            regulated_activity_term_uri = regulated_activity_pc.mechanism["uri"]
        else:
            regulated_activity_term = "GO:0003674"
            regulated_activity_term_uri = None
        connection_clone = pc.clone()
        connection_clone.regulated_activity["term"] = regulated_activity_term
        if connection_clone.regulated_activity["term"] == None or p_connections.contains(connection_clone, check_ref=True):
            skipped_count += 1
            continue
        else:
            pc.regulated_activity["term"] = regulated_activity_term
            pc.regulated_activity["uri"] = regulated_activity_term_uri
            # pc.individuals[pc.regulated_activity["term"]] = regulated_activity_term_uri

        # model = pc.declare_entities(model)

        # enabled_by_stmt_a = model.writer.emit(model.individuals[pc.mechanism_go_term], ENABLED_BY, model.individuals[pc.full_id_a()])
        # if pc.mechanism["term"] in pc.individuals and not model_contains_statement(model, pc.individuals[pc.mechanism["term"]], ENABLED_BY, pc.full_id_a()):
        full_statement = pc.full_statement_bnode_in_model(model)
        # if pc.mechanism["term"] in pc.individuals and full_statement is None:
        if full_statement is None:
            print("Hey " + pc.full_id_a())
            # if pc.id_a == "Q13315" and pc.id_b == "P38398":
            #     print("Dang " + pc.pmid[0])
            model = pc.declare_entities(model)

            enabled_by_stmt_a_triple = (pc.mechanism["uri"], ENABLED_BY, pc.individuals[pc.full_id_a()])
            if enabled_by_stmt_a_triple in model.writer.writer.graph:
                enabled_by_stmt_a = next(model.writer.writer.graph.triples(enabled_by_stmt_a_triple))
            else:
                enabled_by_stmt_a = model.writer.emit(enabled_by_stmt_a_triple[0], enabled_by_stmt_a_triple[1], enabled_by_stmt_a_triple[2])
                axiom_a = model.add_axiom(enabled_by_stmt_a)
            enabled_by_stmt_b_triple = (pc.regulated_activity["uri"], ENABLED_BY, pc.individuals[pc.full_id_b()])
            if enabled_by_stmt_b_triple in model.writer.writer.graph:
                enabled_by_stmt_b = next(model.writer.writer.graph.triples(enabled_by_stmt_b_triple))
            else:
                enabled_by_stmt_b = model.writer.emit(enabled_by_stmt_b_triple[0], enabled_by_stmt_b_triple[1], enabled_by_stmt_b_triple[2])
                axiom_b = model.add_axiom(enabled_by_stmt_b)

            # Connect the two activities
            # source_id = model.individuals[pc.mechanism_go_term]
            try:
                source_id = pc.mechanism["uri"]
            except KeyError as err:
                pc.print()
                print(pc.individuals)
                raise err
            property_id = URIRef(expand_uri(pc.relation))
            # target_id = model.individuals[pc.regulated_activity_term]
            target_id = pc.regulated_activity["uri"]
            if not model_contains_statement(model, source_id, property_id, pc.regulated_activity["term"]):
                # Annotate source MF GO term NamedIndividual with relation code-target MF term URI
                model.writer.emit(source_id, property_id, target_id)
                # Add axiom (Source=MF term URI, Property=relation code, Target=MF term URI)
                relation_axiom = model.writer.emit_axiom(source_id, property_id, target_id)
                model.add_evidence(relation_axiom, "EXP", ["PMID:" + pmid for pmid in pc.pmid])
        else:
            print("2")

        # pc.print()

    with open(model.filepath, 'wb') as f:
        model.writer.writer.serialize(destination=f)

    print(skipped_count)

if __name__ == '__main__':
    main()
    print("hey")