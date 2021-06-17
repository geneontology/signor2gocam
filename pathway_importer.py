from gocamgen.gocamgen import GoCamModel, GoCamEvidence
from ontobio.vocabulary.relations import OboRO
from rdflib.term import URIRef
from rdflib.namespace import Namespace, OWL
from prefixcommons.curie_util import expand_uri
from pathway_connections import PathwayConnectionSet
from entity_models import SignorProtein, SignorMicroRNA, SignorSmallMolecule
from util import OntologyTerm
import argparse
import datetime

ro = OboRO()
ENABLED_BY = URIRef(expand_uri(ro.enabled_by))
HAS_INPUT = URIRef(expand_uri(OntologyTerm.HAS_INPUT.value))
HAS_OUTPUT = URIRef(expand_uri(OntologyTerm.HAS_OUTPUT.value))
EXP_ECO_CODE = "ECO:0000269"

parser = argparse.ArgumentParser()
parser.add_argument('-f', "--filename", type=str, required=True,
                    help="Input filename of SIGNOR pathway data")
parser.add_argument('-t', "--model_title", nargs='+',
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


def pathway_connection_filter_protein_binding(p_connections):
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
                # pc_list.connections.remove(the_good_one)
                # p_connections.remove_list(pc_list.connections)
                p_connections.remove_connection(the_bad_one)
    return p_connections


def generate_model(filename, title):
    model = GoCamModel(title)

    p_connections = PathwayConnectionSet.parse_file(filename)
    linenum = 1

    total_pcs = len(p_connections.connections)
    print(total_pcs, "initial pathway_connections")

    p_connections = pathway_connection_filter_protein_binding(p_connections)

    # fill in regulated activities
    for pc in p_connections.connections:
        # Setup
        pc.declare_a_to_mechanism(model, EXP_ECO_CODE)

    # Now that the a's are declared, go check on the b's.
    for pc in p_connections.connections:
        # Look for triples "anything" -enabled_by-> entity B
        entity_b_pcs = p_connections.find_by_id_a(pc.id_b())
        # If doesn't exist, declare entity B and "anything" becomes root MF, then emit enabled_by
        # TODO
        # Emit reg relation from mechanism URI to entity B triples' activities
        mechanism_uri = pc.mechanism["uri"]
        regulatory_relation = pc.relation
        evidence = pc.gocam_evidence(EXP_ECO_CODE)
        if len(entity_b_pcs) == 0:
            # BPC was likely filtered out due to BPC.entity B not being acceptable type (e.g. phenotype)
            # Declare pc.entity B? A and B should be valid by this point
            # A protein -> B complex
            # A complex -> B phenotype
            # Create and load new PC
            print("No downstream pathway_connections for", pc)
        for bpc in entity_b_pcs:
            participant_relation = HAS_INPUT
            is_small_mol_catalysis = False
            # catalytic activity
            if pc.mechanism["term"] == "GO:0003824" and pc.b_is_small_mol():
                is_small_mol_catalysis = True
                if pc.effect.startswith("up-regulates"):
                    participant_relation = HAS_OUTPUT
            if not pc.a_is_small_mol():
                # mechanism -has_input/output-> entity_b
                has_input_triple = (mechanism_uri, participant_relation, bpc.entity_a.uri)
                if len(model.triples_by_ids(*has_input_triple)) == 0:
                    model.writer.emit(*has_input_triple)
                has_input_axiom = model.find_or_create_axiom(*has_input_triple)
                model.add_evidence(has_input_axiom, evidence=evidence)
                if is_small_mol_catalysis:
                    # Skip adding causal relation
                    continue

            # Add intermediary biological process for these regulatory mechanisms
            intermediary_bp = None
            intermediary_relation = None
            downstream_relation = None
            # ubiquitin protein ligase activity
            if pc.mechanism["term"] == "GO:0061630" and pc.effect.startswith("down-regulates"):
                intermediary_bp = "GO:0043161"  # proteasome-mediated ubiquitin-dependent protein catabolic process
                intermediary_relation = OntologyTerm.POSITIVELY_REGULATES
                downstream_relation = OntologyTerm.NEGATIVELY_REGULATES
            # transcription regulator activity
            if pc.mechanism["term"] == "GO:0140110":
                intermediary_bp = "GO:0009299"  # mRNA transcription
                if pc.effect.startswith("down-regulates"):
                    intermediary_relation = OntologyTerm.NEGATIVELY_REGULATES
                else:
                    intermediary_relation = OntologyTerm.POSITIVELY_REGULATES
                downstream_relation = OntologyTerm.POSITIVELY_REGULATES
            # mRNA 3'-UTR binding
            if pc.mechanism["term"] == "GO:0003730":
                if isinstance(pc.entity_a, SignorMicroRNA) and pc.effect.startswith("down-regulates"):
                    intermediary_bp = "GO:0035195"  # gene silencing by miRNA
                    intermediary_relation = OntologyTerm.POSITIVELY_REGULATES
                    downstream_relation = OntologyTerm.NEGATIVELY_REGULATES
                elif isinstance(pc.entity_a, SignorProtein):
                    intermediary_bp = "GO:0000956"  # nuclear-transcribed mRNA catabolic process
                    if pc.effect.startswith("down-regulates"):
                        intermediary_relation = OntologyTerm.NEGATIVELY_REGULATES
                    else:
                        intermediary_relation = OntologyTerm.POSITIVELY_REGULATES
                    downstream_relation = OntologyTerm.NEGATIVELY_REGULATES
            if intermediary_bp:
                # Extend the statement a bit
                intermediary_bp_uri = model.declare_individual(intermediary_bp)
                # mechanism -has_input-> entity_b
                has_input_triple = (intermediary_bp_uri, HAS_INPUT, bpc.entity_a.uri)
                model.writer.emit(*has_input_triple)
                has_input_axiom = model.writer.emit_axiom(*has_input_triple)
                model.add_evidence(has_input_axiom, evidence=evidence)
                # downstream relation (intermediary_bp -?-> regulated_activity) is static for some of these
                intermediary_triple = (mechanism_uri, URIRef(expand_uri(intermediary_relation.value)), intermediary_bp_uri)
                model.writer.emit(*intermediary_triple)
                intermediary_axiom = model.writer.emit_axiom(*intermediary_triple)
                model.add_evidence(intermediary_axiom, evidence=evidence)
                mechanism_uri, regulatory_relation = intermediary_bp_uri, downstream_relation

            # mechanism -regulates-> regulated_activity OR mechanism -regulates-> intermediary BP -regulates-> regulated_activity
            regulated_activity_uri = bpc.mechanism["uri"]
            regulation_triple = (mechanism_uri, URIRef(expand_uri(regulatory_relation.value)), regulated_activity_uri)
            model.writer.emit(*regulation_triple)
            regulation_axiom = model.writer.emit_axiom(*regulation_triple)
            model.add_evidence(regulation_axiom, evidence=evidence)

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
        model_title = " ".join(args.model_title)
    else:
        model_title = args.outfile
    
    model = generate_model(args.filename, model_title)
    model.write(args.outfile)

if __name__ == '__main__':
    main()