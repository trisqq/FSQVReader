from src.base import SequenceReader


class FastaReader(SequenceReader):
    """Читатель файлов в формате FASTA

    Наследует базовый интерфейс SequenceReader и реализует методы
    проверки, чтения и анализа последовательностей FASTA
    """
    def __init__(self,filename):
        super().__init__(filename) #self.__filename=filename

    def is_fasta(self):
        """Проверяет, соответствует ли файл формату FASTA

        Файл считается корректным FASTA, если каждая последовательность
        начинается со строки-заголовка '>' и содержит хотя бы одну
        непустую строку с данными до следующего заголовка

        Возвращает:
            bool: True - если файл валиден как FASTA, иначе - False
                  Также возвращает False при ошибке чтения файла (OSError)
        """

        try:
            with open(self._filename, "r", encoding="utf-8") as file:
                has_header = False
                has_sequence = False

                for line in file:
                    line = line.strip()

                    if not line:
                        continue

                    if line.startswith(">"):
                        if has_header and not has_sequence:
                            return False

                        has_header = True
                        has_sequence = False

                    else:
                        if not has_header:
                            return False

                        has_sequence = True

                return has_header and has_sequence

        except OSError:
            return False


    def validate_sequence(self, sequence):
        """Проверяет, что последовательность состоит только из нуклеотидов

        Допустимые символы: A, T, U, G, C.
        Пробелы и переводы строк игнорируются

        Возвращает:
            bool: True, если все символы допустимы, если нет -  False
        """
        symbols = sequence.replace(" ", "").replace("\n", "").upper()
        nucleotide = set("ATUGC")
        protein = set("ACDEFGHIKLMNPQRSTVWY")

        for symbol in symbols:
            if (symbol not in nucleotide) and (symbol not in protein):
                return False
        return True

    def read(self):
        """Читает FASTA-файл и выдаёт записи по одной

        Многострочные последовательности склеиваются в одну строку
        ID — первое слово после '>'

        Возвращает:
            запись с полями id (str) и data (str)
        """
        with open(self._filename, 'r', encoding="utf-8") as file:
            seq_id = None
            seq_data = []

            for line in file:
                is_header, text = self._parse_line(line)

                if is_header:
                    if seq_id is not None:
                        yield {"id": seq_id, "data": "".join(seq_data)}

                    seq_id = text
                    seq_data = []

                else:
                    if text:
                        seq_data.append(text)

            if seq_id is not None:
                yield {"id": seq_id, "data": "".join(seq_data)}

    def get_sequence(self, seq_id):
        """Возвращает последовательность по идентификатору
            str/None: последовательность или None, если ID не найден
        """
        for record in self.read():
            if record["id"] == seq_id:
                return record["data"]
        return None #Информация не найдена

    def get_sequence_count(self):
        """Возвращает количество записей в файле
            int: число последовательностей (0 для пустого файла)
        """
        count = 0
        for i in self.read():
            count += 1
        return count

    def get_average_length(self):
        """Вычисляет среднюю длину последовательностей в файле

        Возвращает:
            float: средняя длина записи, 0 - если записей нет
        """
        whole_length = 0
        count = 0
        for record in self.read():
            whole_length += len(record["data"]) #функция считает длину строки
            count += 1

        if count > 0:
            return whole_length / count
        else:
            return 0

    def read_sequences(self):
        """Возвращает множество идентификаторов записей
            set[str]: множество ID (по факту - из одного элемента)
        """
        result = {}
        for record in self.read():
            result[record["id"]] = record["data"]
        return result

    def _parse_line(self, line):
        """Разбирает строку FASTA-файла

        Возвращает:
            tuple[bool, str]: (True, id) - для заголовка '>',
            (False, "") - для пустой строки,
            (False, данные) - для строки с последовательностью
        """
        line = line.strip()

        if not line:
            return False, ""

        if line.startswith(">"):
            return True, line[1:].split()[0]

        return False, line
