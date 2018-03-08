import csv
from gocamgen.gocamgen import GoCamModel
from signor_complex import SignorComplexFactory
from ontobio.vocabulary.relations import OboRO
from rdflib.term import URIRef
from prefixcommons.curie_util import expand_uri

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

class PathwayConnection():
    def __init__(self, id_a, id_b, mechanism, effect, direct, relation, pmid, linenum):
        self.id_a = id_a
        self.id_b = id_b
        self.effect = effect
        self.direct = direct
        self.relation = relation
        self.pmid = pmid
        self.linenum = linenum

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
        # Classes
        if self.full_id_a() not in model.classes:
            model.declare_class(self.full_id_a())
        if self.full_id_b() not in model.classes:
            model.declare_class(self.full_id_b())

        # Individuals
        # if self.full_id_a() not in model.individuals:
        if self.full_id_a() not in self.individuals:
            uri_a = model.declare_individual(self.full_id_a())
            # model.individuals[self.full_id_a()] = uri_a
            self.individuals[self.full_id_a()] = uri_a
        # if self.full_id_b() not in model.individuals:
        if self.full_id_b() not in self.individuals:
            uri_b = model.declare_individual(self.full_id_b())
            # model.individuals[self.full_id_b()] = uri_b
            self.individuals[self.full_id_b()] = uri_b

        self.mechanism["uri"] = model.declare_individual(self.mechanism["term"])
        self.individuals[self.mechanism["term"]] = self.mechanism["uri"]
        if self.regulated_activity["uri"] is None:
            self.regulated_activity["uri"] = model.declare_individual(regulated_activity_term)
        self.individuals[self.regulated_activity["term"]] = self.regulated_activity["uri"]


        return model

    def full_id_a(self):
        return self.full_id(self.id_a)
    def full_id_b(self):
        return self.full_id(self.id_b)

    def full_id(self, id):
        return "UniProtKB:" + id

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

        a_enables_triples = []
        for id_a in model.uri_list_for_individual(self.full_id_a()):
            for mech_uri in model.uri_list_for_individual(self.mechanism["term"]):
                if (mech_uri, ENABLED_BY, id_a) in graph:
                    a_enables_triples.append((mech_uri, ENABLED_BY, id_a))

        b_enables_triples = []
        for id_b in model.uri_list_for_individual(self.full_id_b()):
            for reg_act in model.uri_list_for_individual(self.regulated_activity["term"]):
                if (reg_act, ENABLED_BY, id_b) in graph:
                    b_enables_triples.append((reg_act, ENABLED_BY, id_b))

        for a_triple in a_enables_triples:
            for b_triple in b_enables_triples:
                candidate_reg_triple = (a_triple[0], URIRef(expand_uri(self.relation)), b_triple[0])
                if candidate_reg_triple in graph:
                    return candidate_reg_triple

class PathwayConnectionSet():
    def __init__(self):
        self.connections = []

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

def find_by_id_a(pc_list, id):
    for pc in pc_list:
        if pc.id_a == id:
            return pc

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

model = GoCamModel("test.ttl")
p_connections = PathwayConnectionSet()
linenum = 1
complex_csv_filename = "SIGNOR_complexes.csv"
complexes = SignorComplexFactory(complex_csv_filename).complexes

with open("SIGNOR-G2-M_trans_02_03_18.tsv", "r") as f:
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
        
        # pc.individuals[pc.mechanism["term"]] = model.declare_individual(pc.mechanism["term"]) ### TODO Why do I have to declare this here?
        if not (pc.id_a.startswith("SIGNOR") or pc.id_b.startswith("SIGNOR")):
            p_connections.append(pc)

total_pcs = len(p_connections.connections)
print(total_pcs)
skipped_count = 0

# fill in regulated activities
for pc in p_connections.connections:
    if pc.id_a.startswith("SIGNOR") or pc.id_b.startswith("SIGNOR"):
        # for now to see how model first looks - skip complexes
        continue
    regulated_activity_pc = find_by_id_a(p_connections.connections, pc.id_b)
    if regulated_activity_pc is not None:
        regulated_activity_term = regulated_activity_pc.mechanism["term"]
        # regulated_activity_term_uri = regulated_activity_pc.individuals[regulated_activity_pc.mechanism["term"]]
        regulated_activity_term_uri = regulated_activity_pc.mechanism["uri"]
    else:
        regulated_activity_term = "GO:0003674"
        # regulated_activity_term_uri = model.declare_individual(regulated_activity_term) ### TODO Why do I have to declare this here?
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
        enabled_by_stmt_a = model.writer.emit(pc.individuals[pc.mechanism["term"]], ENABLED_BY, pc.individuals[pc.full_id_a()])
        axiom_a = model.add_axiom(enabled_by_stmt_a)
        # enabled_by_stmt_b = model.writer.emit(model.individuals[pc.regulated_activity_term], ENABLED_BY, model.individuals[pc.full_id_b()])
        enabled_by_stmt_b = model.writer.emit(pc.individuals[pc.regulated_activity["term"]], ENABLED_BY, pc.individuals[pc.full_id_b()])
        axiom_b = model.add_axiom(enabled_by_stmt_b)
        

        # Connect the two activities
        # source_id = model.individuals[pc.mechanism_go_term]
        try:
            source_id = pc.individuals[pc.mechanism["term"]]
        except KeyError as err:
            pc.print()
            print(pc.individuals)
            raise err
        property_id = URIRef(expand_uri(pc.relation))
        # target_id = model.individuals[pc.regulated_activity_term]
        target_id = pc.individuals[pc.regulated_activity["term"]]
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