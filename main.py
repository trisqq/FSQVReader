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
    
    
if __name__ == "__main__":
    main()