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
EXP_ECO_CODE = "ECO:0000269"

parser = argparse.ArgumentParser()
parser.add_argument('-f', "--filename", type=str, required=True,
                    help="Input filename of SIGNOR pathway data")
parser.add_argument('-t', "--model_title",
                    help="Model title. Defaults to --outfile value.")
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

def generate_model(filename, title):
    model = GoCamModel(title)

    p_connections = PathwayConnectionSet.parse_file(filename)
    linenum = 1

    total_pcs = len(p_connections.connections)
    print(total_pcs, "initial pathway_connections")
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

    # fill in regulated activities
    for pc in p_connections.connections:

        # Declare entity A and mechanism
        pc.declare_a(model)
        pc.mechanism["uri"] = model.declare_individual(pc.mechanism["term"])
        # Emit mechanism -enabled_by -> entity_a
        pc.enabled_by_stmt_a = model.writer.emit(pc.mechanism["uri"], ENABLED_BY, pc.entity_a.uri)
        contributors = []
        if pc.annotator:
            contributors = [pc.annotator]
        evidence = GoCamEvidence(EXP_ECO_CODE, ["PMID:" + pmid for pmid in pc.references],
                                 contributors=contributors)
        axiom_a = model.add_axiom(pc.enabled_by_stmt_a, evidence=evidence)
        # model.add_evidence(axiom_a, "EXP", ["PMID:" + pmid for pmid in pc.references])

    # Now that the a's are declared, go check on the b's.
    for pc in p_connections.connections:
        # Look for triples "anything" -enabled_by-> entity B
        entity_b_pcs = p_connections.find_by_id_a(pc.id_b())
        # If doesn't exist, declare entity B and "anything" becomes root MF, then emit enabled_by
        # TODO
        # Emit reg relation from mechanism URI to entity B triples' activities
        mechanism_uri = pc.mechanism["uri"]
        if len(entity_b_pcs) == 0:
            print("No downstream pathway_connections for", pc)
        for bpc in entity_b_pcs:
            # mechanism -has_input-> entity_b
            has_input_triple = (mechanism_uri, HAS_INPUT, bpc.entity_a.uri)
            model.writer.emit(*has_input_triple)
            has_input_axiom = model.writer.emit_axiom(*has_input_triple)

            # mechanism -regulates-> regulated_activity
            regulated_activity_uri = bpc.mechanism["uri"]
            regulation_triple = (mechanism_uri, URIRef(expand_uri(pc.relation)), regulated_activity_uri)
            model.writer.emit(*regulation_triple)
            regulation_axiom = model.writer.emit_axiom(*regulation_triple)
            if pc.annotator:
                contributors = [pc.annotator]
            evidence = GoCamEvidence(EXP_ECO_CODE, ["PMID:" + pmid for pmid in pc.references],
                                     contributors=contributors)
            model.add_evidence(regulation_axiom, evidence=evidence)

    print(skipped_count, "causal statements skipped")
    print(len(p_connections.connections), "pathway_connections at finish")

    grouped = map(lambda x:x.id_a, p_connections.connections)
    print(grouped)

    return model
    

def main():

    ## Organize connection objects
    ## Declare entity A GPs and MFs
    ## Add "Has_input" relations between MF and entity B GPs
    ##      If entity B not declared, declare it
    ## Connect regulation relations to all MF's enabled by entity B
    ##      If no MF for entity B, add root MF enabled by B

    args = parser.parse_args()

    if args.model_title:
        model_title = args.model_title
    else:
        model_title = args.outfile
    
    model = generate_model(args.filename, model_title)
    model.write(args.outfile)

if __name__ == '__main__':
    main()