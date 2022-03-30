
class House():
    def __init__(self, square_footage, price):
        self.square_footage = square_footage
        self.price = price

    def cost_per_square_foot(self):
        return self.square_footage / self.price


from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

record = SeqRecord(Seq("MKQHKAMIVALIVIVITAVVAALVTRK", IUPAC.protein),
                   id="YP_02529.1", name="HokC", description="toxic membrane protein, small")

my_dna = Seq("“AGTACA”")
my_dna.reverse_complement()