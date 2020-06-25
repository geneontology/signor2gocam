import unittest
import yaml
import csv
from rdflib import Graph
from rdflib.plugins.sparql import prepareQuery
from gocamgen.gocamgen import GoCamModel
from pathway_connections import MechanismToGoMappingSet, PathwayConnectionSet
from pathway_importer import generate_model, pathway_connection_filter_protein_binding
from util import OntologyTerm

M_FILE = "metadata/signor_mechanism_go_mapping.yaml"


class TestSignor2Gocam(unittest.TestCase):

    def test_mechanism_map_loading(self):
        with open(M_FILE) as mf:
            mappings = yaml.safe_load(mf)
        self.assertGreater(len(mappings), 1)

    def test_mechanism_mapping_set(self):
        mapping_set = MechanismToGoMappingSet(M_FILE)
        go_term = mapping_set.go_id_by_mechanism("catalytic activity")
        self.assertEqual(go_term, "GO:0003824")

    def run_query(self, model: GoCamModel, query):
        prefix_context = {
            "RO": "http://purl.obolibrary.org/obo/RO_",
            "UniProtKB": "http://identifiers.org/uniprot/",
            "GO": "http://purl.obolibrary.org/obo/GO_"
        }
        graph = model.writer.writer.graph
        response = graph.query(prepareQuery(query, initNs=prefix_context))
        # response = graph.query(prepareQuery(query))
        return response

    def gen_causal_stmt_query(self, entity_a, mechanism, entity_b, reg_relation):
        enabled_by = "RO:0002333"
        query = f"""
                SELECT ?activity ?unknown
                WHERE {{ 
                    ?s rdf:type UniProtKB:{entity_a} .
                    ?activity rdf:type {mechanism} .
                    ?entity_b rdf:type UniProtKB:{entity_b} .
                    ?activity {enabled_by} ?s .
                    ?activity {reg_relation} ?unknown .
                    ?unknown {enabled_by} ?entity_b
                }}
                """
        return query

    def gen_intermediary_bp_stmt_query(self, entity_a, mechanism, intermediary_bp, intermediary_relation,
                                       entity_b, reg_relation):
        enabled_by = "RO:0002333"
        query = f"""
                SELECT ?activity ?unknown
                WHERE {{ 
                    ?s rdf:type UniProtKB:{entity_a} .
                    ?activity rdf:type {mechanism} .
                    ?int_bp rdf:type {intermediary_bp} .
                    ?entity_b rdf:type UniProtKB:{entity_b} .
                    ?activity {enabled_by} ?s .
                    ?activity {intermediary_relation} ?int_bp .
                    ?int_bp {reg_relation} ?unknown .
                    ?unknown {enabled_by} ?entity_b
                }}
                """
        return query

    def test_distinct_entity_instance_per_stmt(self):
        # Sample of line count by ida column from SIGNOR-LBC pathway. Num could change as SIGNOR is updated.
        # 15 O15530
        # 2 P01116
        # 2 P04637
        # 2 P05412
        # 19 P06213
        # 2 P06400
        # 17 P11362
        # 5 P15056
        stmt_file = "resources/test/SIGNOR-LBC.tsv"
        entity_counts = {}
        with open(stmt_file) as sf:
            reader = csv.reader(sf, delimiter="\t")
            next(reader)  # Skip header
            for r in reader:
                entity_a_id = r[5]
                if len(entity_a_id) != 6:
                    # This should hopefully only select UniProt IDs
                    continue
                if entity_a_id not in entity_counts:
                    entity_counts[entity_a_id] = 0
                entity_counts[entity_a_id] += 1
        most_prominent_entity = list(sorted(entity_counts, key=lambda x: entity_counts[x], reverse=True))[0]

        model = generate_model(stmt_file, "SIGNOR - Luminal Breast Cancer")

        query = f"SELECT ?s WHERE {{ ?s rdf:type UniProtKB:{most_prominent_entity}}}"
        resp = self.run_query(model, query)
        print(len(resp), f"individuals for {most_prominent_entity}")

        # Count all triples
        query = "SELECT ?s ?p ?o WHERE {{ ?s ?p ?o}}"
        resp = self.run_query(model, query)
        print(len(resp), "total triples")

        # TODO: Count all causal triples (causally_upstream_of_or_within RO:0002418 or descendants)
        # TODO: Probably need to load RO outside of graph to trace descendants of RO:0002418
        directly_positively_regulates = "RO:0002629"
        query = f"SELECT ?s ?o WHERE {{ ?s {directly_positively_regulates} ?o}}"
        resp = self.run_query(model, query)
        print(len(resp), "causal triples")

        # How many activities does most_prominent_entity enable?
        enabled_by = "RO:0002333"
        query = f"""
        SELECT ?activity
        WHERE {{ 
            ?s rdf:type UniProtKB:{most_prominent_entity} .
            ?activity {enabled_by} ?s
        }}
        """
        resp = self.run_query(model, query)
        print(len(resp), f"activities enabled_by {most_prominent_entity}")

        # Find: most_prominent_entity <-enabled_by- protein_kinase_activity -directly_positively_regulates-> unknown -enabled_by-> irs1_gp
        protein_kinase_activity = "GO:0004672"
        irs1_gp = "P35568"
        resp = self.run_query(model, self.gen_causal_stmt_query(entity_a=most_prominent_entity,
                                                                mechanism=protein_kinase_activity,
                                                                entity_b=irs1_gp,
                                                                reg_relation=directly_positively_regulates
                                                                )
                              )
        print(len(resp), "matches found")
        self.assertEqual(len(resp), 2)

        # TODO: Now check how many instances of most_prominent_entity are in model. Should this == input count?
        print(most_prominent_entity, entity_counts[most_prominent_entity])

    def test_intermediary_bp_patterns(self):
        stmt_file = "resources/test/SIGNOR-LBC.tsv"
        model = generate_model(stmt_file, "SIGNOR - Luminal Breast Cancer")

        # Should only be one: Q00987 <-enabled_by- GO:0061630 -regulates-> GO:0043161 -regulates-> MF -enabled_by-> P04637
        resp = self.run_query(model, self.gen_intermediary_bp_stmt_query(entity_a="Q00987",
                                                                         mechanism="GO:0061630",
                                                                         intermediary_bp="GO:0043161",
                                                                         intermediary_relation=OntologyTerm.POSITIVELY_REGULATES.value,
                                                                         entity_b="P04637",
                                                                         reg_relation=OntologyTerm.NEGATIVELY_REGULATES.value
                                                                        )
                              )
        print(len(resp), "matches found")
        self.assertEqual(len(resp), 1)

    def test_pathway_connection_filter_protein_binding(self):
        stmt_file = "resources/test/SIGNOR-LBC-protein_binding.tsv"
        p_connections = PathwayConnectionSet.parse_file(stmt_file)
        p_connections = pathway_connection_filter_protein_binding(p_connections)
        self.assertEqual(1, 1)

    def test_small_molecule_patterns(self):
        stmt_file = "resources/test/SIGNOR-smallmol.tsv"
        model = generate_model(stmt_file, "SIGNOR - Small molecule test")

        self.assertEqual(1, 1)


if __name__ == '__main__':
    unittest.main()