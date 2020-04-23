import csv
import itertools
import yaml
from ontobio.vocabulary.relations import OboRO
from rdflib.term import URIRef
from prefixcommons.curie_util import expand_uri
from signor_complex import SignorComplexFactory
from naming_conventions import NamingConvention

ro = OboRO()

ENABLED_BY = URIRef(expand_uri(ro.enabled_by))


class MechanismToGoMapping:
    def __init__(self, mechanism, mi_id, go_id, relation):
        self.mechanism = mechanism
        self.mi_id = mi_id
        self.go_id = go_id
        self.relation = relation


class MechanismToGoMappingSet:
    def __init__(self, mapping_file=None):
        self.mappings = []
        if mapping_file:
            with open(mapping_file) as mf:
                mappings = yaml.load(mf)
            for m in mappings:
                self.mappings.append(MechanismToGoMapping(
                    mechanism=m["MECHANISM"],
                    mi_id=m["MI_ID"],
                    go_id=m["GO_ID"],
                    relation=m["RELATION"]
                ))

    def go_id_by_mechanism(self, mechanism):
        for m in self.mappings:
            if m.mechanism == mechanism:
                return m.go_id


class SignorEntity:
    pass


class SignorProtein(SignorEntity):
    pass


# TODO: Refactor `PathwayConnection`
# * Connect causal statements together in networkx graph
# 	* This will reduce need to query RDF triples
# * Then write out to rdflib
class PathwayConnection:
    MECHANISM_GO_MAPPING = MechanismToGoMappingSet("metadata/signor_mechanism_go_mapping.yaml")
    complex_csv_filename = "SIGNOR_complexes.csv"
    COMPLEXES = SignorComplexFactory(complex_csv_filename).complexes

    def __init__(self, id_a, id_b, mechanism, effect, direct: bool, relation, references: list, linenum=None):
        self.id_a = id_a  # TODO: Turn these into entity fields for SignorProtein or SignorComplex objects
        self.id_b = id_b
        self.effect = effect
        self.direct = direct
        self.relation = relation
        self.references = references
        self.linenum = linenum
        self.complex_a = self.complex_from_id(self.id_a)
        self.complex_b = self.complex_from_id(self.id_b)
        # by default mechanism = molecular function
        mechanism_term = "GO:0003674"
        if self.direct:
            if mechanism:
                mechanism_term = self.MECHANISM_GO_MAPPING.go_id_by_mechanism(mechanism)
        self.mechanism = {
            "name": mechanism,
            "uri": None,
            "term": mechanism_term
        }
        self.regulated_activity = {
            "name": None,
            "uri": None,
            "term": None
        }

        self.individuals = {}
        self.enabled_by_stmt_a = None

    @staticmethod
    def parse_line(line: dict, linenum: int=None):
        direct = False
        if line["DIRECT"] in ["YES", "t"]:
            direct = True
        relation = PathwayConnection.determine_relation(effect=line["EFFECT"], direct=direct)

        pc = PathwayConnection(
            id_a=line["IDA"],
            id_b=line["IDB"],
            mechanism=line["MECHANISM"],
            effect=line["EFFECT"],
            direct=direct,
            relation=relation,
            references=[line["PMID"]],
            linenum=linenum
        )
        return pc

    @classmethod
    def determine_relation(cls, effect: str, direct: bool):
        # If up-regulates (including any variants of this), use RO:0002629 if DIRECT,
        # and use RO:0002213 if not DIRECT or UNKNOWN
        relation = None
        if effect.startswith("up-regulates"):
            if direct:
                relation = "RO:0002629"
            else:
                relation = "RO:0002213"
        # If down-regulates (including any variants of this), use RO:0002630 if DIRECT,
        # and use RO:0002212 if not DIRECT or UNKNOWN
        elif effect.startswith("down-regulates"):
            if direct:
                relation = "RO:0002630"
            else:
                relation = "RO:0002212"
        # If unknown, use RO:0002211
        elif effect in ["unknown", ""]:
            relation = "RO:0002211"

        return relation

    @classmethod
    def complex_from_id(cls, entity_id: str):
        if NamingConvention.is_complex(entity_id) and entity_id in cls.COMPLEXES:
            return cls.COMPLEXES[entity_id]

    def print(self):
        print(f"[UniProtKB:{self.id_a}] <- enabled_by – [{self.mechanism['term']}] – [{self.relation}]-> [{self.regulated_activity['term']}] – enabled_by-> [UniProtKB:{self.id_b}]")

    def declare_entities(self, model):
        self.declare_a(model)
        self.declare_b(model)

        return model

    def declare_a(self, model):
        # Class
        if self.class_id_a() not in model.classes:
            model.declare_class(self.class_id_a())

        # Individuals
        # if self.full_id_a() not in self.individuals:
        if self.full_id_a() not in model.individuals:
            if self.a_is_complex():
                uri_a = self.complex_a.declare_entities(model)
            else:
                uri_a = model.declare_individual(self.full_id_a())
            self.individuals[self.full_id_a()] = uri_a
        else:
            self.individuals[self.full_id_a()] = model.individuals[self.full_id_a()]

        # self.mechanism["uri"] = model.declare_individual(self.mechanism["term"])  # Segregate from singular entity declaration
        # self.individuals[self.mechanism["term"]] = self.mechanism["uri"]

        return model

    def declare_b(self, model):
        # Class
        if self.class_id_b() not in model.classes:
            model.declare_class(self.class_id_b())

        # Individuals
        if self.full_id_b() not in self.individuals and self.regulated_activity["uri"] is None:
            if self.b_is_complex():
                uri_b = self.complex_b.declare_entities(model)
            else:
                uri_b = model.declare_individual(self.full_id_b())
            self.individuals[self.full_id_b()] = uri_b
            self.regulated_activity["uri"] = model.declare_individual(self.regulated_activity["term"])
        else:
            for t in model.writer.writer.graph.triples((self.regulated_activity["uri"],ENABLED_BY,None)):
                self.individuals[self.full_id_b()] = t[2]

        self.individuals[self.regulated_activity["term"]] = self.regulated_activity["uri"]

        return model

    def full_id_a(self):
        return NamingConvention.full_id(self.id_a)
    def full_id_b(self):
        return NamingConvention.full_id(self.id_b)
    def class_id_a(self):
        return NamingConvention.class_id(self.id_a)
    def class_id_b(self):
        return NamingConvention.class_id(self.id_b)

    def a_is_complex(self):
        return NamingConvention.is_complex(self.id_a)
    def b_is_complex(self):
        return NamingConvention.is_complex(self.id_b)

    def clone(self):
        new_connection = PathwayConnection(self.id_a, self.id_b, self.mechanism["name"], self.effect, self.direct, self.relation, self.references, self.linenum)
        new_connection.mechanism = self.mechanism
        return new_connection

    def equals(self, pathway_connection, check_ref=False):
        if self.id_a == pathway_connection.id_a and self.id_b == pathway_connection.id_b and self.mechanism == pathway_connection.mechanism and self.relation == pathway_connection.relation and self.regulated_activity == pathway_connection.regulated_activity:
            if check_ref:
                if set(self.references) == set(pathway_connection.references):
                    return True
            else:
                return True
        return False

    def full_statement_bnode_in_model(self, model):
        # Find all existing URI's for IDA, IDB, mech, and reg. Check if statements exist for these URI combos. Might need SPARQL or further triple querying refinement (e.g. triple annotated with "owl:NamedIndividual")
        # mechanism["term"] ENABLED_BY self.id_a
        # regulated_activity["term"] ENABLED_BY self.id_b
        # mechanism["term"] REGULATES regulated_activity["term"]
        graph = model.writer.writer.graph

        # a_enables_triples = []
        # for id_a in model.uri_list_for_individual(self.full_id_a()):
        #     for mech_uri in model.uri_list_for_individual(self.mechanism["term"]):
        #         if (mech_uri, ENABLED_BY, id_a) in graph:
        #             a_enables_triples.append((mech_uri, ENABLED_BY, id_a))
        a_enables_triples = model.triples_by_ids(self.mechanism["term"], ENABLED_BY, self.full_id_a())

        # b_enables_triples = []
        # for id_b in model.uri_list_for_individual(self.full_id_b()):
        #     for reg_act in model.uri_list_for_individual(self.regulated_activity["term"]):
        #         if (reg_act, ENABLED_BY, id_b) in graph:
        #             b_enables_triples.append((reg_act, ENABLED_BY, id_b))
        b_enables_triples = model.triples_by_ids(self.regulated_activity["term"], ENABLED_BY, self.full_id_b())

        for a_triple in a_enables_triples:
            for b_triple in b_enables_triples:
                candidate_reg_triple = (a_triple[0], URIRef(expand_uri(self.relation)), b_triple[0])
                if candidate_reg_triple in graph:
                    return candidate_reg_triple


def upper_first(iterator):
    return itertools.chain([next(iterator).upper()], iterator)


class PathwayConnectionSet():
    def __init__(self):
        self.connections = []

    @staticmethod
    def parse_file(filename):
        pc_set = PathwayConnectionSet()

        linenum = 0
        if filename:
            with open(filename, "r") as f:
                data = list(csv.DictReader(upper_first(f), delimiter="\t"))
                for line in data:
                    linenum += 1

                    acceptable_types = ['protein', 'complex']
                    if line["TYPEA"] not in acceptable_types or \
                       line["TYPEB"] not in acceptable_types or \
                       line["EFFECT"] == "form complex":
                        continue

                    pc = PathwayConnection.parse_line(line, linenum=linenum)
                    pc_set.add(pc)

        return pc_set


    def add(self, pathway_connection):
        existing_connection = self.find(pathway_connection)
        if existing_connection:
            # Causal statement already exists so just add reference
            existing_connection.references = set(existing_connection.references) | set(pathway_connection.references)
        else:
            self.connections.append(pathway_connection)

    def contains(self, pathway_connection, check_ref=False):
        for connection in self.connections:
            if connection.equals(pathway_connection, check_ref=check_ref):
                return True
        return False

    def find(self, pathway_connection, check_ref=False):
        for connection in self.connections:
            if connection.equals(pathway_connection, check_ref=check_ref):
                return connection

    def find_by_id_a(self, id):
        pcs = []
        for pc in self.connections:
            if pc.id_a == id:
                pcs.append(pc)
        return pcs

    def find_other_regulated_activity(self, id_b):
        regulated_pcs = self.find_by_id_a(id_b)
        filtered_reg_pcs = []
        for pc in regulated_pcs:
            if pc.mechanism["term"] != "GO:0003674":
                filtered_reg_pcs.append(pc)
        if len(filtered_reg_pcs) > 0:
            return filtered_reg_pcs[0]

    def find_all_by_id_a_and_id_b(self, pathway_connection):
        found_connections = PathwayConnectionSet()
        for pc in self.connections:
            if pc.id_a == pathway_connection.id_a and pc.id_b == pathway_connection.id_b:
                found_connections.add(pc)
        return found_connections

    def find_by_mech_term(self, term):
        for pc in self.connections:
            if pc.mechanism["term"] == term:
                return pc

    def remove_list(self, pc_list):
        new_connection_list = self.connections
        for dead_pc in pc_list:
            new_connection_list = [pc for pc in new_connection_list if not pc.equals(dead_pc)]
        self.connections = new_connection_list