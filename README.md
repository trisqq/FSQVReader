# FSQVReader

Читалки для биоинформатических форматов данных: **FASTA**, **FASTQ**, **SAM** и **VCF**.

Библиотека позволяет получать количество записей, среднюю длину последовательностей,
заголовки файлов и их группы, а также статистику по хромосомам и выборки
по заданному геномному отрезку. Для FASTQ строятся графики качества в стиле FastQC.

## Установка и запуск

Нужен Python 3.12 или новее.

```bash
# 1. Склонировать репозиторий
git clone https://github.com/trisqq/FSQVReader.git
cd FSQVReader

# 2. Создать и активировать окружение conda
conda create -n fsqvreader python=3.12 -y
conda activate fsqvreader

# 3. Установить зависимости
pip install -r requirements.txt

# 4. Запустить демонстрацию
python main.py
