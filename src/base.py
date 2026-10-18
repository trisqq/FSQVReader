import abc

class Reader(abc.ABC):
    def __init__ (self, filename):
        self._filename = filename
        self.file = open(filename, mode = "r")

    def read(self):
        """Считывает строчки в файле, выдает их в виде генератора"""
        for line in self.file:
            yield self._parse_line(line)

    def close(self):
        """Закрывает файл"""
        self.file.close()

    @abc.abstractmethod
    def _parse_line(self, line):
        """Определяет тип строки в файле"""
        pass



class GenomicDataReader(Reader):

    @abc.abstractmethod
    def get_chromosomes(self):
        pass

    @abc.abstractmethod
    def get_reference_genomes(self):
        pass

    @abc.abstractmethod
    def validate_coordinates(self, chr, pos):
        pass



