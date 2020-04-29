import csv
import itertools
import yaml
from typing import List
from ontobio.vocabulary.relations import OboRO
from rdflib.term import URIRef
from prefixcommons.curie_util import expand_uri
from entity_factories import SignorEntityFactory
from entity_models import SignorEntity

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
                mappings = yaml.safe_load(mf)
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


class AnnotatorOrcidMapping:
    def __init__(self, annotator_name, orcid: str):
        self.annotator_name = annotator_name
        self.orcid = orcid

    def full_orcid_uri(self):
        orcid_prefix = "http://orcid.org/"
        if self.orcid.startswith(orcid_prefix):
            return self.orcid
        return orcid_prefix + self.orcid

class AnnotatorOrcidMappingSet:
    def __init__(self, file=None):
        self.mappings = []
        if file:
            with open(file) as mf:
                mappings = list(csv.DictReader(mf, delimiter="\t"))
                for m in mappings:
                    self.mappings.append(AnnotatorOrcidMapping(
                        annotator_name=m["ANNOTATOR"],
                        orcid=m["ORCID"],
                    ))

    def orcid_by_name(self, annotator_name):
        for m in self.mappings:
            if m.annotator_name == annotator_name:
                return m.full_orcid_uri()

# TODO: Refactor `PathwayConnection`
# * Connect causal statements together in networkx graph
# 	* This will reduce need to query RDF triples
# * Then write out to rdflib
class PathwayConnection:
    MECHANISM_GO_MAPPING = MechanismToGoMappingSet("metadata/signor_mechanism_go_mapping.yaml")
    ANNOTATOR_ORCID_MAPPING = AnnotatorOrcidMappingSet("metadata/annotator_orcid.tsv")
    def __init__(self, entity_a: SignorEntity, entity_b: SignorEntity, mechanism, effect, direct: bool, relation, references: list, annotator, linenum=None):
        self.entity_a = entity_a
        self.entity_b = entity_b
        self.effect = effect
        self.direct = direct
        self.relation = relation
        self.references = references
        self.linenum = linenum
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

        self.annotator = None
        if annotator:
            self.annotator = self.ANNOTATOR_ORCID_MAPPING.orcid_by_name(annotator)

        self.individuals = {}
        self.enabled_by_stmt_a = None

    @staticmethod
    def parse_line(line: dict, linenum: int=None):
        entity_a = SignorEntityFactory.determine_entity(entity_id=line["IDA"], entity_name=line["ENTITYA"])
        entity_b = SignorEntityFactory.determine_entity(entity_id=line["IDB"], entity_name=line["ENTITYB"])

        direct = False
        if line["DIRECT"] in ["YES", "t"]:
            direct = True
        relation = PathwayConnection.determine_relation(effect=line["EFFECT"], direct=direct)

        pc = PathwayConnection(
            entity_a=entity_a,
            entity_b=entity_b,
            mechanism=line["MECHANISM"],
            effect=line["EFFECT"],
            direct=direct,
            relation=relation,
            references=[line["PMID"]],
            annotator=line["ANNOTATOR"],
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
        # If unknown, use RO:0002211 (regulates)
        elif effect in ["unknown", ""]:
            relation = "RO:0002211"

        return relation

    def __str__(self):
        return f"[UniProtKB:{self.id_a()}] <- enabled_by – [{self.mechanism['term']}] – [{self.relation}]-> [{self.regulated_activity['term']}] – enabled_by-> [UniProtKB:{self.id_b()}]"

    def print(self):
        print(self)

    def declare_entities(self, model):
        self.declare_a(model)
        self.declare_b(model)

    def declare_a(self, model):
        self.entity_a.declare(model)

    def declare_b(self, model):
        self.entity_b.declare(model)

    def id_a(self):
        return self.entity_a.id

    def id_b(self):
        return self.entity_b.id

    @staticmethod
    def _full_id(entity: SignorEntity):
        return entity.full_id()

    def full_id_a(self):
        return self._full_id(self.entity_a)

    def full_id_b(self):
        return self._full_id(self.entity_b)

    @staticmethod
    def _class_id(entity: SignorEntity):
        return entity.class_id()

    def class_id_a(self):
        return self._class_id(self.entity_a)

    def class_id_b(self):
        return self._class_id(self.entity_b)

    @staticmethod
    def _is_complex(entity: SignorEntity):
        return entity.is_complex()

    def a_is_complex(self):
        return self._is_complex(self.entity_a)

    def b_is_complex(self):
        return self._is_complex(self.entity_b)

    def clone(self):
        new_connection = PathwayConnection(self.entity_a, self.entity_b, self.mechanism["name"], self.effect,
                                           self.direct, self.relation, self.references, self.linenum)
        new_connection.mechanism = self.mechanism
        return new_connection

    def equals(self, pathway_connection, check_ref=False):
        if self.entity_a == pathway_connection.entity_a and self.entity_b == pathway_connection.entity_b and \
                self.mechanism == pathway_connection.mechanism and self.relation == pathway_connection.relation \
                and self.regulated_activity == pathway_connection.regulated_activity:
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


class PathwayConnectionSet:
    def __init__(self):
        self.connections: List[PathwayConnection] = []

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

    def add(self, pathway_connection: PathwayConnection):
        existing_connection = self.find(pathway_connection)
        if existing_connection:
            # Causal statement already exists so just add reference
            existing_connection.references = set(existing_connection.references) | set(pathway_connection.references)
        else:
            self.connections.append(pathway_connection)

    def contains(self, pathway_connection: PathwayConnection, check_ref=False):
        for connection in self.connections:
            if connection.equals(pathway_connection, check_ref=check_ref):
                return True
        return False

    def find(self, pathway_connection: PathwayConnection, check_ref=False):
        for connection in self.connections:
            if connection.equals(pathway_connection, check_ref=check_ref):
                return connection

    def find_by_id_a(self, id) -> List[PathwayConnection]:
        pcs = []
        for pc in self.connections:
            if pc.id_a() == id:
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
            if pc.id_a() == pathway_connection.id_a() and pc.id_b() == pathway_connection.id_b():
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
