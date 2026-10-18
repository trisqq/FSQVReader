import abc

class Reader(abc.ABC):
    def __init__ (self, filename):
        self._filename = filename
        self.file = open(filename, mode = "r")

    def read(self):
        for line in self.file:
            yield self._parse_line(line)

    def close(self):
        self.file.close()

    @abc.abstractmethod
    def _parse_line(self, line):
        pass