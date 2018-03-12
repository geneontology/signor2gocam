import csv
from naming_conventions import NamingConvention
from rdflib.term import URIRef, Literal
from rdflib.namespace import RDFS
from prefixcommons.curie_util import expand_uri

complexes = []

# class SignorComplex():
class SignorGrouping():
    def __init__(self, signor_id, name, entities):
        self.id = signor_id
        self.name = name
        self.entities = entities

class SignorComplex(SignorGrouping):
    def declare_entities(self, model):
        uri = model.declare_individual("GO:0032991")
        model.writer.writer.graph.add((uri, RDFS.label, Literal(str(self.name))))
        for entity in self.entities:
            entity_full_id = NamingConvention.full_id(entity)
            entity_uri = model.declare_individual(entity_full_id)
            part_of_stmt = model.writer.emit(uri, URIRef(expand_uri("BFO:0000051")), entity_uri)
            model.add_axiom(part_of_stmt)
            "uri BFO:0000051 entity_uri"
        return uri

class SignorProteinFamily(SignorGrouping):
    def declare_entities(self, model):
        uri = model.declare_individual("GO:0032991")
        model.writer.writer.graph.add((uri, RDFS.label, Literal(str(self.name))))
        for entity in self.entities:
            entity_full_id = NamingConvention.full_id(entity)
            entity_uri = model.declare_individual(entity_full_id)
            part_of_stmt = model.writer.emit(uri, URIRef(expand_uri("BFO:0000051")), entity_uri)
            model.add_axiom(part_of_stmt)
            "uri BFO:0000051 entity_uri"
        return uri

class SignorGroupingFactory():
    NAME_FIELD = None
    GROUPING_CLASS = None

    def __init__(self, filename):
        self.grouping = {}
        with open(filename, "r") as f:
            data = list(csv.DictReader(f, delimiter=";"))


            for line in data:
                entities = []

                for entity in line['LIST OF ENTITIES'].split(", "):
                    entities.append(entity.strip())

                args = {"signor_id" : line['SIGNOR ID'], "name" : line[self.NAME_FIELD], "entities" : entities}
                sig_grouping = eval(self.GROUPING_CLASS)(**args)
                self.grouping[sig_grouping.id] = sig_grouping

class SignorComplexFactory(SignorGroupingFactory):
    def __init__(self, filename):
        self.NAME_FIELD = "COMPLEX NAME"
        self.GROUPING_CLASS = "SignorComplex"
        SignorGroupingFactory.__init__(self, filename)
        self.complexes = self.grouping

class SignorProteinFamilyFactory(SignorGroupingFactory):
    def __init__(self, filename):
        self.NAME_FIELD = "PROT. FAMILY NAME"
        self.GROUPING_CLASS = "SignorProteinFamily"
        SignorGroupingFactory.__init__(self, filename)
        self.families = self.grouping

def main():
    complex_list = SignorComplexFactory("SIGNOR_complexes.csv")
    complex = complex_list.complexes["SIGNOR-C87"]
    complex_cc = "GO:0032991"
    # create individual for go term
    for entity in complex.entities:
        relation = "complex_cc has_part " + entity
        # create individual for entity
        # state relation
        print(relation)
    pathway_line_mech_or_reg_activity = "pc.whatevs"
    # state "pathway_line_mech_or_reg_activity ENABLED_BY complex_cc"

    pf_list = SignorProteinFamilyFactory("SIGNOR_PF.csv")
    pf = pf_list.families["SIGNOR-PF11"]
    print("PF list for " + pf.name + ":")
    for entity in pf.entities:
        print(entity)

if __name__ == "__main__":
    main()