from fastq import FastqReader

if __main__ = "main":
    try:
        fastq = FastqReader("test.fastq")
        print(f"Количество последовательностей: {fastq.get_sequence_count()}")
        print(f"Средняя длина: {fastq.get_average_length():.2f} bp")        
        print("Строим графики (это может занять пару секунд)...")
        fastq.generate_plots(limit=1000) 
            
        fastq.close()
    except FileNotFoundError:
        print("Файл test.fastq не найден. Создайте его для теста.")
