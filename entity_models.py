import abc
from gocamgen.gocamgen import GoCamModel
from naming_conventions import NamingConvention
from rdflib.term import URIRef, Literal
from rdflib.namespace import RDFS
from prefixcommons.curie_util import expand_uri

HAS_PART = URIRef(expand_uri("BFO:0000051"))


class SignorEntity:
    def __init__(self, entity_id, name, uri=None):
        self.id = entity_id
        self.name = name
        self.uri = uri

    def full_id(self):
        return NamingConvention.full_id(self.id)

    def class_id(self):
        return NamingConvention.class_id(self.id)

    def is_complex(self):
        return isinstance(self, SignorComplex)

    def __eq__(self, other):
        if isinstance(other, SignorEntity):
            return self.id == other.id and self.name == other.name
        return False

    def __str__(self):
        return f"{self.id} - {self.name}"

    @abc.abstractmethod
    def declare(self, model: GoCamModel):
        return


class SignorProtein(SignorEntity):
    def declare(self, model):
        self.uri = model.declare_individual(self.full_id())
        return self.uri


class SignorMicroRNA(SignorEntity):
    def declare(self, model):
        self.uri = model.declare_individual(self.full_id())
        return self.uri


class SignorSmallMolecule(SignorEntity):
    def declare(self, model):
        # First, check if instance already exists
        existing_uris = model.uri_list_for_individual(self.full_id())
        if len(existing_uris) > 0:
            self.uri = existing_uris[0]
        else:
            self.uri = model.declare_individual(self.full_id())
        return self.uri

    def full_id(self):
        if self.id.startswith("CHEBI:"):
            return self.id
        else:
            return f"CHEBI:{self.id}"


class SignorGrouping(SignorEntity):
    def __init__(self, signor_id, name, entities):
        SignorEntity.__init__(self, signor_id, name)
        self.entities = entities


class SignorComplex(SignorGrouping):
    def declare(self, model):
        return self.declare_entities(model)

    def declare_entities(self, model):
        self.uri = model.declare_individual("GO:0032991")
        model.writer.writer.graph.add((self.uri, RDFS.label, Literal(str(self.name))))
        for entity in self.entities:
            entity_full_id = NamingConvention.full_id(entity)
            entity_uri = model.declare_individual(entity_full_id)
            part_of_stmt = model.writer.emit(self.uri, HAS_PART, entity_uri)
            model.add_axiom(part_of_stmt)
            "uri BFO:0000051 entity_uri"
        return self.uri

    def uri_in_model(self, model):
        graph = model.writer.writer.graph
        complex_term = "GO:0032991"
        complex_uris = model.uri_list_for_individual(complex_term)
        for c_uri in complex_uris:
            found_entities = [] # Collect connected entities for set comparison
            for entity in self.entities:
                for entity_uri in model.uri_list_for_individual(NamingConvention.full_id(entity)):
                    if (c_uri, HAS_PART, entity_uri) in graph:
                        found_entities.append(entity)
            if set(self.entities) == set(found_entities):
                print(c_uri)
                return c_uri


class SignorProteinFamily(SignorGrouping):
    def declare_entities(self, model):
        uri = model.declare_individual("GO:0032991")
        model.writer.writer.graph.add((uri, RDFS.label, Literal(str(self.name))))
        for entity in self.entities:
            entity_full_id = NamingConvention.full_id(entity)
            entity_uri = model.declare_individual(entity_full_id)
            part_of_stmt = model.writer.emit(uri, HAS_PART, entity_uri)
            model.add_axiom(part_of_stmt)
            "uri BFO:0000051 entity_uri"
        return uri
