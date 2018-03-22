import csv
from ontobio.vocabulary.relations import OboRO
from rdflib.term import URIRef
from prefixcommons.curie_util import expand_uri
from signor_complex import SignorComplexFactory
from naming_conventions import NamingConvention

ro = OboRO()

ENABLED_BY = URIRef(expand_uri(ro.enabled_by))

MECHANISM_GO_MAPPING = {
    "acetylation" : "GO:0061733",
    "binding" : "GO:0005515",
    "chemical activation" : None,
    "chemical inhibition" : None,
    "cleavage" : "GO:0008233",
    "deacetylation" : "GO:0033558",
    "demethylation" : "GO:0032451",
    "dephosphorylation" : "GO:0004721",
    "destabilization" : "GO:0003674",
    "desumoylation" : "GO:0070140",
    "deubiquitination" : "GO:0004843",
    "glycosylation" : "GO:0016757",
    "GAP" : "GO:0005096",
    "GEF" : "GO:0005085",
    "hydroxylation" : "GO:0036140",
    "lipidation" : "GO:0016747",
    "methylation" : "GO:0008276",
    "neddylation" : "GO:0061663",
    "oxidation" : "GO:0003674",
    "palmitoylation" : "GO:0016409",
    "phosphorylation" : "GO:0004672",
    "post transcriptional regulation" : "GO:0035925",
    "relocalization" : "GO:0003674",
    "small molecule catalysis" : None,
    "stabilization" : "GO:0003674",
    "sumoylation" : "GO:0061665",
    "transcriptional activation" : "GO:0003700",
    "transcriptional regulation" : "GO:0003700",
    "transcriptional repression" : "GO:0003700",
    "trimethylation (histone)" : "GO:0003674",
    "tyrosination" : "GO:0004835",
    "ubiquitination" : "GO:0061630"
}

complex_csv_filename = "SIGNOR_complexes.csv"
COMPLEXES = SignorComplexFactory(complex_csv_filename).complexes

class PathwayConnection():
    def __init__(self, id_a, id_b, mechanism, effect, direct, relation, pmid, linenum):
        self.id_a = id_a
        self.id_b = id_b
        self.effect = effect
        self.direct = direct
        self.relation = relation
        self.pmid = pmid
        self.linenum = linenum
        self.complex_a = None
        try:
            if NamingConvention.is_complex(self.id_a) and self.id_a in COMPLEXES:
                self.complex_a = COMPLEXES[self.id_a]
        except TypeError as err:
            print(self.id_a)
            raise err
        self.complex_b = None
        if NamingConvention.is_complex(self.id_b) and self.id_b in COMPLEXES:
            self.complex_b = COMPLEXES[self.id_b]

        if self.direct == "NO":
            mechanism_term = "GO:0003674"
        elif self.direct == "YES":
            mechanism_term = MECHANISM_GO_MAPPING[mechanism]
            # self.mechanism["term"] = MECHANISM_GO_MAPPING["stabilization"]
        self.mechanism = { "name" : mechanism, "uri" : None, "term" : mechanism_term }
        self.regulated_activity = { "name" : None, "uri" : None, "term" : None }

        self.individuals = {}

    def print(self):
        print("[UniProtKB:{ida}] <- enabled_by – [{mechanism}] – [{relation}]-> [{regulated_activity}] – enabled_by-> [UniProtKB:{idb}]".format(ida=self.id_a,
                                                                                                                                                mechanism=self.mechanism["term"],
                                                                                                                                                relation=self.relation,
                                                                                                                                                regulated_activity=self.regulated_activity["term"],
                                                                                                                                                idb=self.id_b))

    def declare_entities(self, model):
        self.declare_a(model)
        self.declare_b(model)

        return model

    def declare_a(self, model):
        # Class
        if self.full_id_a() not in model.classes:
            model.declare_class(self.full_id_a())

        # Individuals
        if self.full_id_a() not in self.individuals:
            if self.a_is_complex():
                uri_a = self.complex_a.declare_entities(model)
            else:
                uri_a = model.declare_individual(self.full_id_a())
            self.individuals[self.full_id_a()] = uri_a

        self.mechanism["uri"] = model.declare_individual(self.mechanism["term"])
        self.individuals[self.mechanism["term"]] = self.mechanism["uri"]

        return model

    def declare_b(self, model):
        # Class
        if self.full_id_b() not in model.classes:
            model.declare_class(self.full_id_b())

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

    def a_is_complex(self):
        return NamingConvention.is_complex(self.id_a)
    def b_is_complex(self):
        return NamingConvention.is_complex(self.id_b)

    def clone(self):
        new_connection = PathwayConnection(self.id_a, self.id_b, self.mechanism["name"], self.effect, self.direct, self.relation, self.pmid, self.linenum)
        new_connection.mechanism = self.mechanism
        return new_connection

    def equals(self, pathway_connection, check_ref=False):
        if self.id_a == pathway_connection.id_a and self.id_b == pathway_connection.id_b and self.mechanism == pathway_connection.mechanism and self.relation == pathway_connection.relation and self.regulated_activity == pathway_connection.regulated_activity:
            if check_ref:
                if set(self.pmid) == set(pathway_connection.pmid):
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

class PathwayConnectionSet():
    def __init__(self, filename):
        self.connections = []
        linenum = 0

        with open(filename, "r") as f:
            data = list(csv.DictReader(f, delimiter="\t"))
            for line in data:
                linenum += 1

                # If up-regulates (including any variants of this), use RO:0002629 if DIRECT, and use RO:0002213 if not DIRECT
                relation = None
                if line["EFFECT"].startswith("up-regulates"):
                    if line["DIRECT"] == "YES":
                        relation = "RO:0002629"
                    elif line["DIRECT"] == "NO":
                        relation = "RO:0002213"
                # If down-regulates (including any variants of this), use RO:0002630 if DIRECT, and use RO:0002212 if not DIRECT
                if line["EFFECT"].startswith("down-regulates"):
                    if line["DIRECT"] == "YES":
                        relation = "RO:0002630"
                    elif line["DIRECT"] == "NO":
                        relation = "RO:0002212"
                # If unknown, use RO:0002211
                if line["EFFECT"] == "unknown":
                    relation = "RO:0002211"
                # If form_complex, ignore these lines for now
                if line["EFFECT"] == "form_complex":
                    continue

                pc = PathwayConnection(
                    line["IDA"],
                    line["IDB"],
                    line["MECHANISM"],
                    line["EFFECT"],
                    line["DIRECT"],
                    relation,
                    [line["PMID"]],
                    linenum
                )

                # if not (pc.id_a.startswith("SIGNOR") or pc.id_b.startswith("SIGNOR") or line["TYPEA"] == "phenotype" or line["TYPEB"] == "phenotype"):
                acceptable_types = ['protein','complex']
                if line["TYPEA"] in acceptable_types and line["TYPEB"] in acceptable_types:
                    if self.find(pc):
                        self.append_reference(pc)
                    else:
                        self.append(pc)


    def append(self, pathway_connection):
        self.connections.append(pathway_connection)

    def append_reference(self, pathway_connection):
        connection = self.find(pathway_connection)
        connection.pmid = set(connection.pmid) | set(pathway_connection.pmid)

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