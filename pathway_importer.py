from gocamgen.gocamgen import GoCamModel, GoCamEvidence
from ontobio.vocabulary.relations import OboRO
from rdflib.term import URIRef
from rdflib.namespace import Namespace, OWL
from prefixcommons.curie_util import expand_uri
from pathway_connections import PathwayConnectionSet
import argparse

ro = OboRO()
ENABLED_BY = URIRef(expand_uri(ro.enabled_by))
HAS_INPUT = URIRef(expand_uri("RO:0002233"))

parser = argparse.ArgumentParser()
parser.add_argument('-f', "--filename", type=str, required=True,
                    help="Input filename of Signor pathway data")
parser.add_argument('-o', "--outfile", type=str, required=True,
                    help="Output filename of generated model")

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

    ## Organize connection objects
    ## Declare entity A GPs and MFs
    ## Add "Has_input" relations between MF and entity B GPs
    ##      If entity B not declared, declare it
    ## Connect regulation relations to all MF's enabled by entity B
    ##      If no MF for entity B, add root MF enabled by B

    args = parser.parse_args()

    model = GoCamModel(args.outfile)
    # p_connections = PathwayConnectionSet("SIGNOR-G2-M_trans_02_03_18.tsv")
    p_connections = PathwayConnectionSet(args.filename)
    linenum = 1
    # complex_csv_filename = "SIGNOR_complexes.csv"
    # complexes = SignorComplexFactory(complex_csv_filename).complexes

    total_pcs = len(p_connections.connections)
    print(total_pcs)
    skipped_count = 0

    # Toss out connections according to precedence rules:
    # protein kinase activity should be chosen over protein binding
    # This should be separate from/before any OWL individuals are declared
    for pc in p_connections.connections:
        # Get list of all pcs with pc.id_a and pc.id_b
        # If len of list is > 1
        #   Look for protein kinase activity (GO:0004672), delete others if remaining are all protein binding (GO:0005515)
        pc_list = p_connections.find_all_by_id_a_and_id_b(pc)
        if len(pc_list.connections) > 1:
            the_good_one = pc_list.find_by_mech_term("GO:0004672")
            the_bad_one = pc_list.find_by_mech_term("GO:0005515")
            if the_good_one is not None and the_bad_one is not None:
                pc_list.connections.remove(the_good_one)
                p_connections.remove_list(pc_list.connections)
            elif the_good_one is not None:
                [uncertain_pc.print() for uncertain_pc in pc_list.connections]
                p_connections.remove_list(pc_list.connections)

    # fill in regulated activities
    for pc in p_connections.connections:
        regulated_activity_pc = p_connections.find_other_regulated_activity(pc.id_b)
        if regulated_activity_pc is not None:
            regulated_activity_term = regulated_activity_pc.mechanism["term"]
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

        # enabled_by_stmt_a_triple = (pc.mechanism["uri"], ENABLED_BY, pc.individuals[pc.full_id_a()])
        if pc.a_is_complex():
            entity_a = pc.complex_a.uri_in_model(model)
            if entity_a is None:
                model = pc.declare_a(model)
                entity_a = pc.complex_a.uri_in_model(model)
        else:
            entity_a = pc.full_id_a()
        # Don't care about existing "statements", just look for existing entity A GP and always create new enabled by statement
        # if entity_a is not None:
        #     enabled_by_stmt_a_triples = model.triples_by_ids(pc.mechanism["term"], ENABLED_BY, entity_a)
        # else:
        #     enabled_by_stmt_a_triples = []
        # if len(enabled_by_stmt_a_triples) == 0:

        # If triple A doesn't exist for entities, declare individuals and create it
        # enabled_by_stmt_a = model.writer.emit(pc.mechanism["term"], ENABLED_BY, pc.full_id_a())
        if entity_a is None or not isinstance(entity_a, URIRef):
            for uri_a in model.uri_list_for_individual(pc.full_id_a()):
                pc.individuals[pc.full_id_a()] = uri_a
            if pc.full_id_a() not in pc.individuals:
                model = pc.declare_a(model)
        else:
            pc.individuals[pc.full_id_a()] = entity_a
        pc.mechanism["uri"] = model.declare_individual(pc.mechanism["term"])
        pc.individuals[pc.mechanism["term"]] = pc.mechanism["uri"]
        pc.enabled_by_stmt_a = model.writer.emit(pc.mechanism["uri"], ENABLED_BY, pc.individuals[pc.full_id_a()])
        axiom_a = model.add_axiom(pc.enabled_by_stmt_a, GoCamEvidence("EXP", ["PMID:" + pmid for pmid in pc.pmid]))
        # model.add_evidence(axiom_a, "EXP", ["PMID:" + pmid for pmid in pc.pmid])

        # else:
        #     pc.enabled_by_stmt_a = enabled_by_stmt_a_triples[0]

    # Now that the a's are declared, go check on the b's.
    for pc in p_connections.connections:
        if pc.b_is_complex():
            entity_b = pc.complex_b.uri_in_model(model)
            entity_b_uris = [entity_b]
        else:
            entity_b = pc.full_id_b()
            entity_b_uris = model.uri_list_for_individual(entity_b)
        if entity_b is not None:
            # enabled_by_stmt_b_triples = model.triples_by_ids(pc.regulated_activity["term"], ENABLED_BY, entity_b)
            enabled_by_stmt_b_triples = model.triples_by_ids(None, ENABLED_BY, entity_b)
        else:
            enabled_by_stmt_b_triples = []
        # If pointing to activity
        regulated_activity_uris = []
        for b_triple in enabled_by_stmt_b_triples:
            regulated_activity_uris.append(b_triple[0])

        if len(regulated_activity_uris) == 0:
            model = pc.declare_b(model)
            enabled_by_stmt_b = model.writer.emit(pc.regulated_activity["uri"], ENABLED_BY,
                                                  pc.individuals[pc.full_id_b()])
            axiom_b = model.add_axiom(enabled_by_stmt_b)
            # model.add_evidence(axiom_b, "EXP", ["PMID:" + pmid for pmid in pc.pmid])  # Maybe don't want to add evidence to B since we're assuming these statements?
            regulated_activity_uris.append(enabled_by_stmt_b[0])
            # entity_b_uris.append(pc.individuals[pc.full_id_b()])
            entity_b_uris = [pc.individuals[pc.full_id_b()]]

        for entity_b_uri in entity_b_uris:
            model.writer.emit(pc.enabled_by_stmt_a[0], HAS_INPUT, entity_b_uri)
            relation_axiom = model.writer.emit_axiom(pc.enabled_by_stmt_a[0], HAS_INPUT, entity_b_uri)

        # Connect the two activities
        # Decouple this from ENABLED_BY statements to allow multiple regulation relations from one GP-MF node - issue #2
        source_id = pc.enabled_by_stmt_a[0]
        property_id = URIRef(expand_uri(pc.relation))
        # if not model_contains_statement(model, source_id, property_id, pc.regulated_activity["term"]):  # Make into for loop
        for reg_activity_uri in regulated_activity_uris:
            target_id = reg_activity_uri
            # Annotate source MF GO term NamedIndividual with relation code-target MF term URI
            model.writer.emit(source_id, property_id, target_id)
            # Add axiom (Source=MF term URI, Property=relation code, Target=MF term URI)
            relation_axiom = model.writer.emit_axiom(source_id, property_id, target_id)
            model.add_evidence(relation_axiom, "EXP", ["PMID:" + pmid for pmid in pc.pmid])

    model.write(args.outfile)

    print(skipped_count)

    grouped = map(lambda x:x.id_a, p_connections.connections)
    print(grouped)

if __name__ == '__main__':
    main()