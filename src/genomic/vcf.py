import re

import pandas as pd

from ..base import GenomicDataReader


class Variant:
    def __init__(self, chrom, pos, id, ref, alt, qual, filter, info, samples):
        self.chrom = chrom
        self.pos = pos
        self.id = id
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter = filter
        self.info = info
        self.samples = samples

    def __repr__(self):
        return f"Variant: {self.__dict__}"


class Genotype:
    def __init__(self, sample, fields):
        self.sample = sample
        self.fields = fields

    def __repr__(self):
        return f"Genotype: {self.__dict__}"


class VcfHeader:
    def __init__(self, lines):
        self.lines = lines

    def get_group(self, group):
        """Возвращает записи одной группы заголовка (INFO, FILTER, FORMAT, contig и т.д.) в виде словарей"""
        groups = []
        for line in self.lines:
            if line.startswith(f"##{group}="):
                content = line.split("=", 1)[1]
                groups.append(dict(re.findall(r'(\w+)=("[^"]*"|[^,>]+)', content)))
        return groups


class VcfReader(GenomicDataReader):

    def __init__(self, filename):
        super().__init__(filename)
        self.samples = []
        self.filename = filename

    def _parse_line(self, line):
        line = line.rstrip('\n')
        if not line:
            return None
        if line.startswith('##'):
            return line
        if line.startswith('#CHROM'):
            self.samples = line.split('\t')[9:]
            return line
        chrom, pos, vid, ref, alt, qual, flt, info, fmt, *sample_cols = line.split('\t')
        keys = fmt.split(':')
        samples = {name: dict(zip(keys, col.split(':'))) for name, col in zip(self.samples, sample_cols)}
        qual_value = None if qual == '.' else float(qual)
        return Variant(chrom, int(pos), vid, ref, alt, qual_value, flt, info, samples)

    def read_variants(self):
        """Возвращает генератор вариантов, все записи в память не сохраняются"""
        self.file.seek(0)
        for record in self.read():
            if isinstance(record, Variant):
                yield record

    def get_header(self):
        """Возвращает заголовок VCF файла"""
        self.file.seek(0)
        lines = []
        for line in self.file:
            if not line.startswith('#'):
                break
            lines.append(line.rstrip('\n'))
        return VcfHeader(lines)

    def get_variant_count(self):
        """Возвращает общее количество вариантов"""
        return sum(1 for _ in self.read_variants())

    def get_chromosome_stats(self):
        """Возвращает DataFrame: количество вариантов на каждой хромосоме"""
        counts = {}
        for variant in self.read_variants():
            counts[variant.chrom] = counts.get(variant.chrom, 0) + 1
        return pd.DataFrame(list(counts.items()), columns=['chrom', 'count'])

    def get_chromosomes(self):
        return [contig.get('ID') for contig in self.get_header().get_group('contig')]

    def get_reference_genomes(self):
        for line in self.get_header().lines:
            if line.startswith('##reference='):
                return line.split('=', 1)[1]
        return None

    def validate_coordinates(self, chr, pos):
        for contig in self.get_header().get_group('contig'):
            if contig.get('ID') == chr:
                return 0 < pos <= int(contig.get('length', 0))
        return False

    def get_variants_in_region(self, chrom, start, end):
        """Возвращает варианты внутри отрезка chrom:start-end (границы включительно, как в bedtools intersect)"""
        for variant in self.read_variants():
            if variant.chrom == chrom and start <= variant.pos <= end:
                yield variant

    def filter_by_quality(self, min_qual: float):
        """Возвращает варианты с QUAL не ниже min_qual (варианты без QUAL, то есть '.', пропускаются)"""
        for variant in self.read_variants():
            if variant.qual is not None and variant.qual >= min_qual:
                yield variant

    def get_genotype(self, sample: str, variant: Variant):
        """Возвращает генотип образца sample для варианта variant"""
        if sample not in variant.samples:
            return None
        return Genotype(sample, variant.samples[sample])
