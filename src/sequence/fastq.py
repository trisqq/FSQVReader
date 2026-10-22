from src.base import Reader
import os
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

class FastqReader(Reader):
    def _parse_line(self, line):
        return line.strip()

    def read(self):

        """
        Генератор, который читает файл FASTQ блоками по 4 строки.
        """

        while True:
            header = self.file.readline().strip()
            if not header:
                break
            seq = self.file.readline().strip()
            plus = self.file.readline().strip()
            qual = self.file.readline().strip()
            yield 
            {
                'id': header[1:].split()[0], 
                'description': header,
                'sequence': seq,
                'quality': qual
            }

    def get_sequence_count(self):

        """
        Получение количества последовательностей
        """

        self.file.seek(0)
        return sum(1 for _ in self.read())

    def get_average_length(self):

        """
        Получение средней длины последовательностей
        """

        self.file.seek(0)
        total_len = 0
        count = 0
        for record in self.read():
            total_len += len(record['sequence'])
            count += 1
        return total_len / count if count > 0 else 0

    def get_quality_scores(self, seq_id: str) -> list[int]:

        """
        Возвращает список числовых значений качества для заданного ID
        """

        self.file.seek(0)
        for record in self.read():
            if record['id'] == seq_id:
                return [ord(char) - 33 for char in record['quality']]
        return []

    def get_average_quality(self, seq_id: str) -> float:

        """
        Возвращает среднее качество прочтения для заданного ID
        """

        scores = self.get_quality_scores(seq_id)
        return sum(scores) / len(scores) if scores else 0.0
    
    def generate_plots(self, limit=10000, output_dir="plots"):

        """
        Генерирует графики качества и состава.
        """

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        self.file.seek(0)
        per_base_quality = defaultdict(list)
        per_base_content = defaultdict(lambda: {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'N': 0})
        sequence_lengths = []
        count = 0
        for record in self.read():
            if count >= limit:
                break
            seq = record['sequence']
            qual = record['quality']
            sequence_lengths.append(len(seq))
            for i, (nuc, q_char) in enumerate(zip(seq, qual)):
                per_base_quality[i].append(ord(q_char) - 33) # Phred+33
                if nuc.upper() in per_base_content[i]:
                    per_base_content[i][nuc.upper()] += 1
                else:
                    per_base_content[i]['N'] += 1
            count += 1
        self._plot_quality(per_base_quality, output_dir)
        self._plot_content(per_base_content, output_dir)
        self._plot_length_distribution(sequence_lengths, output_dir)
        
        print(f"Графики сохранены в папку {output_dir}/")

    def _plot_quality(self, per_base_quality, output_dir):

        """
        График sample_sam: Per base sequence quality
        """

        positions = sorted(per_base_quality.keys())
        avg_qualities = [sum(per_base_quality[p])/len(per_base_quality[p]) for p in positions]
        plt.figure(figsize=(10, 5))
        sns.lineplot(x=positions, y=avg_qualities, color="blue")
        plt.title("Per base sequence quality")
        plt.xlabel("Position in read (bp)")
        plt.ylabel("Average Quality Score (Phred)")
        plt.axhline(y=20, color='r', linestyle='--', label='Low Quality Threshold (20)')
        plt.axhline(y=30, color='g', linestyle='--', label='Good Quality Threshold (30)')
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "per_base_quality.png"))
        plt.close()

    def _plot_content(self, per_base_content, output_dir):

        """
        График 2: Per base sequence content
        """

        positions = sorted(per_base_content.keys())
        bases = ['A', 'T', 'G', 'C', 'N']
        colors = {'A': 'green', 'T': 'red', 'G': 'blue', 'C': 'orange', 'N': 'gray'}
        plt.figure(figsize=(10, 5))
        for base in bases:
            percentages = []
            for p in positions:
                total = sum(per_base_content[p].values())
                percentages.append((per_base_content[p][base] / total) * 100 if total > 0 else 0)
            sns.lineplot(x=positions, y=percentages, label=base, color=colors[base])
        plt.title("Per base sequence content")
        plt.xlabel("Position in read (bp)")
        plt.ylabel("Percentage (%)")
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "per_base_content.png"))
        plt.close()

    def _plot_length_distribution(self, lengths, output_dir):

        """
        График 3: Sequence length distribution
        """

        plt.figure(figsize=(10, 5))
        sns.histplot(lengths, bins=50, kde=True, color="purple")
        plt.title("Sequence length distribution")
        plt.xlabel("Read length (bp)")
        plt.ylabel("Count")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "length_distribution.png"))
        plt.close()
