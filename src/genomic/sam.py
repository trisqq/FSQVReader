from ..base import GenomicDataReader


class Alignment:
    def __init__(self):
        pass
    def __repr__(self):
        return f"Alignment: {self.__dict__}"



class SamReader(GenomicDataReader):

    def __init__(self, filename):
        super().__init__(filename)
        self.header = []
        self.filename = filename

    def _parse_line(self, line):
        if line.startswith('@'):
            self.header.append(line)
            return self.header
        else:
            marks = ["qname", "flag", "rname", "pos", "mapq", "cigar", "rnext", "pnext", "tlen", "seq", "qual"]
            record = dict(zip(marks, line.split("\t")))
            return record


    def get_chromosomes(self):
        chromosomes = set()
        for y in self.read_alignments():
            chromosomes.add(y.rname)
        return chromosomes


    def get_reference_genomes(self):
        pass

    def validate_coordinates(self, chr, pos):
        result = False
        chromosome = {}
        for y in self.get_chr():
            if chr in y:
                result = True
                chromosome = y
        if  result == True and 0 < pos <= int(chromosome[chr]):
            result = True
        else:
            result = False
        return result



    def read_alignments(self):
        """Возвращает список всех выравниваний"""
        al_list = []
        for e in self.read():
            if isinstance(e, dict):
                alignm = Alignment()
                for k, v in e.items():
                    setattr(alignm, k, v)
                al_list.append(alignm)
        return al_list






    def flag_analyze(self,flag):
        """Анализирует флаг выравнивания, возвращает список с битами"""
        bit_template = [2048, 1024, 512, 256, 128, 64, 32, 16, 8, 4, 2, 1]
        bit_list = []
        for i in bit_template:
            check = flag - i
            if check >= 0:
                bit_list.append(i)
                flag = check
        return bit_list


    def get_header(self):
        """выводит заголовок"""
        new_header = []
        for e in self.read():
            if isinstance(e, list):
                new_header = e

        return new_header


    def filter_alignments(self, filter: int):
        """фильтрует выравнивания по конкретному биту"""
        filtered_alignments = []
        for x in self.read_alignments():
            if filter in self.flag_analyze(int(x.flag)):
                filtered_alignments.append(x)
        return filtered_alignments


    def calculate_coverage(self):
        pass


    def get_chr(self):
        chr_list = []
        for i in self.get_header():
            chr_name = ""
            chr_len = ""
            if i.startswith("@SQ"):
                chr_data = {}
                parts = i.strip().split("\t")
                for part in parts:

                    if part.startswith("SN:"):
                        chr_name = part[3:]
                    elif part.startswith("LN:"):
                        chr_len = part[3:]
                chr_data.update({chr_name: chr_len})
                chr_list.append(chr_data)
        return chr_list








a = SamReader(r"C:\Users\bogda\PycharmProjects\PythonProject\FSQVReader\1")
print(a.validate_coordinates("chr1", 200))
