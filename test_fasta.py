import inspect
import pytest

from src.sequence.fasta import FastaReader


#ошибки при проблемах с файлом
EXPECTED_ERRORS = (OSError, ValueError)

# тестовые данные, длины 8, 4 и 6, среднее 6.0
FASTA_TEXT = (
    ">seq1\n"
    "ATGC\n"
    "GGCC\n"
    ">seq2\n"
    "ATGC\n"
    ">seq3\n"
    "AUGCUA\n"
)


@pytest.fixture
def fasta_file(tmp_path):
    """ Создаёт временный FASTA-файл через встроенную 
    функцию tmp_path в изолированную временную папку, 
    записывает туда тестовый текст и возвращает путь к нему """
    path = tmp_path / "sample.fasta"
    path.write_text(FASTA_TEXT, encoding="utf-8")
    return str(path)


@pytest.fixture
def reader(fasta_file):
    """ Берется путь, который вернула функция fasta_file и создается
    готовый объект FastaReader для тестового файла """
    return FastaReader(fasta_file)



  #Определение формата В блоке определения формата я проверяю метод is_fasta на двух сценариях. Первый тест берет наш подготовленный правильный FASTA-файл через фикстуру reader и проверяет, что метод возвращает True. Второй тест проверяет пограничный случай: я создаю через tmp_path обычный текстовый файл без заголовков >, передаю его в новый объект FastaReader и проверяю, что метод вернет False, то есть не пропустит невалидный файл

def test_is_fasta_true(reader):
    """ Файл с заголовками '>' распознаётся как FASTA """
    assert reader.is_fasta() is True



def test_plain_text_is_not_fasta(tmp_path):
    """ Обычный текст без '>' не считается FASTA """
    path = tmp_path / "notes.txt"
    path.write_text("просто текст\nбез заголовка\n", encoding="utf-8")
    assert FastaReader(str(path)).is_fasta() is False



    #Проверка алфавита

def test_validate_nucleotide(reader):
    """ Символы - нуклеотиды ДНК/РНК допустимы """
    assert reader.validate_sequence("ATGCU") is True


def test_validate_protein(reader):
    """ Символы аминокислот допустимы """
    assert reader.validate_sequence("ACDEFGHIKLMNPQRSTVWY") is True


def test_validate_garbage(reader):
    """ Посторонние символы ZX не проходят проверку """
    assert reader.validate_sequence("ATGCZX") is False



    #Чтение записей

def test_read_returns_all_records(reader):
    """ read() выдаёт все 3 записи файла """
    assert len(list(reader.read())) == 3


def test_read_is_generator(reader):
    """ read() - генератор (для больших файлов) """
    assert inspect.isgenerator(reader.read())


def test_read_sequences(reader):
    """ read_sequences() возвращает идентификаторы по порядку """
    assert list(reader.read_sequences()) == ["seq1", "seq2", "seq3"]


def test_get_sequence(reader):
    """ Последовательность находится с помощью идентификатора """
    assert reader.get_sequence("seq2") == "ATGC"


def test_get_sequence_multiline(reader):
    """ Строки одной записи склеиваются в одну последовательность """
    assert reader.get_sequence("seq1") == "ATGCGGCC"



    #Статистика

def test_sequence_count(reader):
    """ Количество последовательностей равно 3 """
    assert reader.get_sequence_count() == 3


def test_average_length(reader):
    """ Средняя длина: (8 + 4 + 6)/3 = 6 """
    assert reader.get_average_length() == pytest.approx(6.0)



    #Обработка ошибок

def test_missing_file(tmp_path):
    """ Несуществующий файл даёт ошибку, а не молчит """
    with pytest.raises(EXPECTED_ERRORS):
        list(FastaReader(str(tmp_path / "no_file.fasta")).read())


def test_directory_instead_of_file(tmp_path):
    """ Папка вместо файла тоже даёт ошибку """
    with pytest.raises(EXPECTED_ERRORS):
        list(FastaReader(str(tmp_path)).read())
