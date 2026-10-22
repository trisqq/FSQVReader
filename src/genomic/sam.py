from FSQVReader.src.base import GenomicDataReader
import pandas as pd

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
            try:
                chromosomes.add(y.rname)
            except AttributeError:
                continue
        return chromosomes


    def get_reference_genomes(self):
        pass

    def validate_coordinates(self, chr, pos):
        result = False
        chrom_len = 0

        ref_data = self.get_ref_chr_data()

        for y in ref_data:
            if chr == y:
                result = True
                chrom_len = int(ref_data[y])
                break

        if result == True and 0 < pos <= chrom_len:
            return True
        else:
            return False


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
            try:
                if filter in self.flag_analyze(int(x.flag)):
                    filtered_alignments.append(x)
            except AttributeError:
                continue
        return filtered_alignments


    def alignments_coordinates(self):
        coords_list = []
        for x in self.read_alignments():
            try:
                start = int(x.pos)
                end = len(x.seq)
            except AttributeError:
                continue
            coordinates = list(range(start, start + end ))
            coords_list.append(coordinates)
        return coords_list



    def calculate_coverage(self, chrom):
        cov_data = {}
        for key, value in self.get_ref_chr_data().items():
            if chrom == key:
                for y in range(1, int(value)+1):
                    cov_count = 0
                    for j in self.alignments_coordinates():
                        if y in j:
                            cov_count += 1
                    cov_data.update({y: cov_count})
        return cov_data



    def get_ref_chr_data(self):
        chr_list = {}
        for i in self.get_header():
            chr_name = ""
            chr_len = ""
            if i.startswith("@SQ"):
                parts = i.strip().split("\t")
                for part in parts:

                    if part.startswith("SN:"):
                        chr_name = part[3:]
                    elif part.startswith("LN:"):
                        chr_len = part[3:]
                chr_list.update({chr_name: chr_len})
        return chr_list

    def get_header_by_tag(self, tag):
        """
        Возвращает список строк заголовка, начинающихся с указанного тега.
        """
        if not tag.startswith('@'):
            tag = '@' + tag

        header_lines = self.get_header()
        result = []

        for line in header_lines:
            if line.startswith(tag):
                result.append(line)

        return result

    def get_alignments_count(self):
        """
        Возвращает общее количество выравниваний в файле.
        """
        self.file.seek(0)
        alignments = self.read_alignments()
        return len(alignments)

    def get_alignment_stats(self):
        self.file.seek(0)
        alignments = self.read_alignments()

        chroms = []
        for align in alignments:
            if hasattr(align, 'rname'):
                chroms.append(align.rname)

        series = pd.Series(chroms)
        counts = series.value_counts()

        stats = pd.DataFrame({'chromosome': counts.index, 'count': counts.values})
        return stats

    def get_alignments_in_region(self, chrom, start, end):
        """
        Возвращает список выравниваний, пересекающихся с регионом [start, end] на хромосоме chrom.
        Координаты start и end включительно (sample_sam-based).
        """
        self.file.seek(0)
        alignments = self.read_alignments()

        result = []

        for align in alignments:
            if getattr(align, 'rname', None) != chrom:
                continue

            try:
                align_start = int(align.pos)
                align_end = align_start + len(align.seq) - 1

            except (AttributeError, ValueError):
                continue

            if align_start <= end and align_end >= start:
                result.append(align)

        return result






