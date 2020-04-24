
class NamingConvention():
    @staticmethod
    def full_id(id):
        if NamingConvention.is_complex(id):
            return id
        return "UniProtKB:" + id

    @staticmethod
    def class_id(id):
        if NamingConvention.is_complex(id):
            return "GO:0032991"  # protein-containing complex
        return "UniProtKB:" + id

    @staticmethod
    def is_complex(id):
        if id.startswith("SIGNOR-C"):
            return True
        else:
            return False

    @staticmethod
    def is_family(id):
        if id.startswith("SIGNOR-PF"):
            return True
        else:
            return False
