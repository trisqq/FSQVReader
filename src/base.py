from abc import ABC, abstractmethod


class Reader(ABC):
    def __init__ (self, filename):
        self._filename = filename
        self.file = open(filename, mode = "r")  # noqa: SIM115

    def read(self):
        """Считывает строчки в файле, выдает их в виде генератора"""
        for line in self.file:
            yield self._parse_line(line)

    def close(self):
        """Закрывает файл"""
        self.file.close()

    @abstractmethod
    def _parse_line(self, line):
        """Определяет тип строки в файле"""
        raise NotImplementedError



class GenomicDataReader(Reader):

    @abstractmethod
    def get_chromosomes(self):
        raise NotImplementedError

    @abstractmethod
    def get_reference_genomes(self):
        raise NotImplementedError

    @abstractmethod
    def validate_coordinates(self, chr, pos):
        raise NotImplementedError
