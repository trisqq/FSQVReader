from src.genomic.vcf import VcfReader
from src.sequence.fasta import FastaReader


def main():
    filename = "sequence.fasta"
    reader = FastaReader(filename)

    is_fasta = reader.is_fasta()
    print(f"Файл является FASTA файлом: {is_fasta}")

    validate_nucleotide = reader.validate_sequence("ATGCU")
    print(f"Валидная ДНК/РНК: {validate_nucleotide}")

    validate_protein = reader.validate_sequence("ACDEFGHIKLMNPQRSTVWY")
    print(f"Валидный белок: {validate_protein}")

    valid_garbage = reader.validate_sequence("ATGCZX")
    print(f"Строка 'ATGCZX' (Мусор) валидна?: {valid_garbage}")

    records = list(reader.read())
    print(f"Записи: {records}")

    sequence = reader.get_sequence("seq2")
    print(f"Последовательность по идентификаору: {sequence}")

    count = reader.get_sequence_count()
    print(f"Количество последовательностей: {count}")

    average_length = reader.get_average_length()
    print(f'Средняя длина последовательности: {average_length}')

    read_sequences = reader.read_sequences()
    print(f"Идентификаторы последовательностей: {read_sequences}")

def demo_vcf(filename):
    reader = VcfReader(filename)
    print(f"INFO заголовок: {reader.get_header().get_group('INFO')}")
    print(f"Хромосомы (contig): {reader.get_chromosomes()}")
    print(f"Количество вариантов: {reader.get_variant_count()}")
    print(f"Статистика по хромосомам:\n{reader.get_chromosome_stats()}")
    print(f"Варианты в chr1:1-500: {list(reader.get_variants_in_region('chr1', 1, 500))}")
    print(f"Варианты с QUAL >= 30: {list(reader.filter_by_quality(30))}")
    reader.close()

if __name__ == "__main__":
    main()
    demo_vcf("sample.vcf")
