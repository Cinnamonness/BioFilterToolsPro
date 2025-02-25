# BioFilterToolsPro

-[russian version](#russian-version)

-[english version](#english-version)

### russian version

**BioFilterToolsPro** — это утилита, предназначенная для работы с последовательностями ДНК, РНК и белков, а также для фильтрации последовательностей FASTQ-файла на основе GC-состава, длины рида и порогового значения среднего качества рида (шкала phred33). 

**bio_files_processor** - дополнительная утилита для работы с некоторыми распространенными форматами биологических данных (fasta-файлы, выходные файлы программы BLAST в формате .txt, геномные аннотации .gbk).

Авторы:
* **Программное обеспечение:** *Карицкая Полина* <cinnamonness@gmail.com>, <br/>
Институт Биоинформатики, Санкт-Петербург, Россия. 
* **Идея, руководитель:** *Никита Ваулин*, *Антон Сидорин* <br/>
Институт Биоинформатики, Санкт-Петербург, Россия.

---
## Оглавление
- [Установка](#установка)
- [Возможности](#возможности)
    - [BioFilterToolsPro](#Возможности-BioFilterToolsPro)
    - [bio_files_processor](#Возможности-bio_files_processor)
- [Системные требования](#системные-требования)
- [Пример использования](#Пример-использования)

---
## Установка

- Используя Git, клонируйте репозиторий на ваш локальный компьютер.

```bash
git clone git@github.com:Cinnamonness/BioFilterToolsPro.git
cd BioFilterToolsPro
```

---

## Возможности

### Возможности BioFilterToolsPro
Утилита представляет две основные функции: 

1. **Операции над последовательностями ДНК и РНК:**
    - Определение типа молекулы (ДНК или РНК)
    - Транскрипция
    - Обратный порядок ("реверсирование")
    - Определение комплементарной последовательности
    - Определение обратной комплементарной последовательности

2. **Фильтрация последовательностей из FASTQ-файла на основе:**
    - GC-состава
    - Длины последовательности
    - Порогового значения среднего качества рида (шкала phred33)

### Возможности bio_files_processor

Утилита представляет три основные функции: 

1. **Работа с fasta-файлами**
    - Конвертация многострочного fasta-файла в однострочный
    - Получение валидного fasta-файла

2. **Работа с выходными файлами BLAST**
    - Запись белков с наилучшим совпадением с базой из выходного файла BLAST
    - Получение .txt файла со списком белков, отсортированных в алфавитном порядке

3. **Работа с геномными аннотациями в формате .gbk**
    - Запись определенного количества генов до и после переданных генов интереса вместе с их белковыми последовательностями (translation)
    - Получение валидного fasta-файла

---

## Системные требования
- Python 3.x
- Необходимые библиотеки из requirements.txt

---

## Пример использования: 

### Пример работы с ДНК последовательностями:
```Python
dna = DNASequence("ATGCGA")
print(dna)
print(dna[1:4])
print(dna[1::])
print(dna.complement())
print(dna.reverse())
print(dna.reverse_complement())
print(dna.transcribe())
```
```Python
Sequence: ATGCGA
Sequence: TGC
Sequence: TGCGA
Sequence: TACGCT
Sequence: AGCGTA
Sequence: TCGCAT
Sequence: AUGCGA
```

### Пример работы с РНК последовательностями:
```Python
rna = RNASequence("AUGCGA")
print(str(rna))        
print(rna[1:4])
print(rna[1::])
print(rna.complement())
print(rna.reverse())
print(rna.reverse_complement())
```
```Python
Sequence: AUGCGA
Sequence: UGC
Sequence: UGCGA
Sequence: UACGCU
Sequence: AGCGUA
Sequence: UCGCAU
```

### Пример работы с белковыми последовательностями:
```Python
aa_seq = AminoAcidSequence("GALNQRHKTTYCC")
print(aa_seq)
print(aa_seq[1:4])
print(aa_seq[1::])
print(aa_seq.classify_aminoacids())
```
```Python
Sequence: GALNQRHKTTYCC
Sequence: ALN
Sequence: ALNQRHKTTYCC
non-polar: Count = 3, Percentage = 23.08%
polar uncharged: Count = 7, Percentage = 53.85%
polar negatively charged: Count = 0, Percentage = 0.00%
polar positively charged: Count = 3, Percentage = 23.08%
```

### Пример работы `filter_fastq`
```Python
if __name__ == "__main__":
    arguments = ("../data/example.fastq", 
             output_fastq='../data/filtered.fastq',
             gc_bounds=(0, 80),
             length_bounds=(0, 500),
             quality_threshold=40)
    result = filter_fastq(*arguments)
```
В результате работы функции `filter_fastq` в директории ./data сохранится файл filtered.fastq с отфильтрованными последовательностями из изначального файла example.fastq. 

***Исходный example_fastq.fastq***
```Python
@K00271:89:HHWWNBBXX:2:1101:23277:1068 1:N:0:CAGATC
NATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAATCTCAGACAACAAATCACAGAGTTAAGTCAGTTTACCGCACAAACTNACCAAGTCGGCGAAACAGAAGGTGGCGAC
+
#AAAFFJJFJJFJJJJF-FFFJJFJJJJFJJJJFJJFAFJF-FFJFJFJJJJ7A<FFJJJFAJF<<-JJJFJJF----FA--7-A-7------7A-7----7FJ----7---7-7A-AF-#A<---7---7A-F)-7AA<--77-)--)-
@K00271:89:HHWWNBBXX:2:1101:16792:1121 1:N:0:CAGATC
NGGAAACCGTCGGTTCTGGTTGTGGAGGCGGTTGGTGGTGGCTGTGGATTGGGAAGATGGTTTGGTAAGTTTTGGGCCGGTTTCAAGAAACTAAGCTGGGCTGGGCTGGGCTGGGCTGGGCTAAGCTGGGCTCACGAACCAAGTAAAGTT
+
#AAFFJJJJJJJJJJJJJJFFJFJJJJJFJJ-FFJFFJFJJJJJJJJFFFJJJ<<JF<FJ<FJJJFFJJJJFJJJJJJJJ-AJJFJJJAJJFJJJJFJJJJJJFFJJJJJJFJJJJJJJJJ-AFJFJJJJJFFAJJ-AFJJJFF-7F<77
@K00271:89:HHWWNBBXX:2:1101:18609:1173 1:N:0:CAGATC
NTGAACCTCCAAATTGATTAGAGAGTCACATGCAAATCACTGTTAAAGTCCGACAAACAGTTTACAAGTAACTTTAATTGAATGCCCGAAAATCTTTTAAATGTTCAGACAATACGATGACGATGACAGTTATAAGCGAAACCACTGAGA
+
#AAFFJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ<FJJJFJFJJJJJJJFJJJJFJJFJJJF
@K00271:89:HHWWNBBXX:2:1101:19898:1226 1:N:0:CAGATC
NTGATGACTGTTGCCAAACAATTTGGGAATTCTAGATGGGATTCGAGTTTAGTTTTGGAGTGAGCCTTATAATTTTGGTTCATCAAGGTCAATAAGGATACACTCCCACATTGGTGTTCATTGGGTTAATTTTGGAGTGCCACTCACACC
+
#AAFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ<FFJJFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJAJFJJJJJJJJJJJFJJFJJJJJJJJJJJJJJJJJJ
```
***Отфильтрованный filtered_sequences.fastq***
```Python
@K00271:89:HHWWNBBXX:2:1101:18609:1173 1:N:0:CAGATC
NTGAACCTCCAAATTGATTAGAGAGTCACATGCAAATCACTGTTAAAGTCCGACAAACAGTTTACAAGTAACTTTAATTGAATGCCCGAAAATCTTTTAAATGTTCAGACAATACGATGACGATGACAGTTATAAGCGAAACCACTGAGA
+
#AAFFJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ<FJJJFJFJJJJJJJFJJJJFJJFJJJF
@K00271:89:HHWWNBBXX:2:1101:19898:1226 1:N:0:CAGATC
NTGATGACTGTTGCCAAACAATTTGGGAATTCTAGATGGGATTCGAGTTTAGTTTTGGAGTGAGCCTTATAATTTTGGTTCATCAAGGTCAATAAGGATACACTCCCACATTGGTGTTCATTGGGTTAATTTTGGAGTGCCACTCACACC
+
#AAFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ<FFJJFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJAJFJJJJJJJJJJJFJJFJJJJJJJJJJJJJJJJJJ
```
### Пример работы `convert_multiline_fasta_to_oneline` из ***bio_files_processor***
```Python
if __name__ == "__main__":
    arguments = ('C:/Users/User/Downloads/example_multiline_fasta.fasta', '')
    result = convert_multiline_fasta_to_oneline(*arguments)
```
***Исходный многострочный example_multiline_fasta.fasta***
```Python
>5S_rRNA::NODE_272_length_223_cov_0.720238:18-129(+)
ACGGCCATAGGACTTTGAAAGCACCGCATCCCGTCCGATCTGCGAAGTTAACCAAGATGCCGCCTGGTTAGTACCATGGTGGGGGACCACATGGGAATCCCT
GGTGCTGTG
>16S_rRNA::NODE_4_length_428221_cov_75.638017:281055-282593(-)
TTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAACAGCTTGCTGTTTCGCTGACGAGTGG
CGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACC
TTCGGGCCTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTTGTTGGTGGGGTAACGGCTCACCAAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCA
GCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGA
AGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTC
CGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCG
GGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAAT
ACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGA
TGTCGACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATT
GACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACGGAAGTTTTCAGAGATGAGAAT
GTGCCTTCGGGAACCGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGT
TGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACA
CACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCA
TGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAA
GAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGGGAACCTGCGGTTGGATC
ACCTCCTT
```
***Полученный однострочный output_fasta.fasta***
```Python
>5S_rRNA::NODE_272_length_223_cov_0.720238:18-129(+)
ACGGCCATAGGACTTTGAAAGCACCGCATCCCGTCCGATCTGCGAAGTTAACCAAGATGCCGCCTGGTTAGTACCATGGTGGGGGACCACATGGGAATCCCTGGTGCTGTG
>16S_rRNA::NODE_4_length_428221_cov_75.638017:281055-282593(-)
TTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAACAGCTTGCTGTTTCGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTTGTTGGTGGGGTAACGGCTCACCAAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACGGAAGTTTTCAGAGATGAGAATGTGCCTTCGGGAACCGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGGGAACCTGCGGTTGGATCACCTCCTT
```

### Пример работы `parse_blast_output` из ***bio_files_processor***
```Python
if __name__ == "__main__":
    arguments = ('C:/Users/User/Downloads/example_blast_results.txt', '')
    result = parse_blast_output(*arguments)
```
С результами работы функции parse_blast_output можно ознакомиться в директории examples. example_blast_results.txt - это исходный файл, который нужно было парсить. parse_blast.txt - файл с отобранными белками. 

### Пример работы `select_genes_from_gbk_to_fasta` из ***bio_files_processor***
```Python
if __name__ == "__main__":
    arguments = ('C:/Users/User/Downloads/example_gbk.gbk', 'sucA', 'btuD_1')
    result = select_genes_from_gbk_to_fasta(*arguments, n_before=2, n_after=2, output_fasta='')
```
С результами работы функции select_genes_from_gbk_to_fasta можно ознакомиться в директории examples. example_gbk.gbk - это исходный файл с геномной аннотацией ***E. coli***. gbk.fasta - fasta-файл с выделенным количеством генов до и после каждого из гена интереса и сохраненной белковой последовательностью (translation). 

---
# BioFilterToolsPro

### english version

**BioFilterToolsPro** — is a utility designed to work with DNA, RNA and protein sequences, as well as to filter sequences in FASTQ file based on GC composition, reed length and the threshold value of the average reed quality (phred33 scale).

**bio_files_processor** - is an additional utility for working with some common biological data formats (FASTA files, output files from the BLAST program in .txt format, genome annotations in .gbk format).

Authors:
* **Software:** *Karitskaya Polina* <cinnamonness@gmail.com>, <br/>
Institute of Bioinformatics, Saint-Petersburg, Russia. 
* **Idea, supervisor:** *Nikita Vaulin*, *Anton Sidorin* <br/>
Institute of Bioinformatics, Saint-Petersburg, Russia.

---
## Table of contents
- [How to install](#how-to-install)
- [Features](#features)
    -[Features of BioFilterToolsPro](#Features-of-BioFilterToolsPro)
    -[Features of bio_files_processor](#Features-of-bio_files_processor)
- [Requirements](#requirements)
- [Example](#Example)

---

## How to install

- Clone the repository to your local machine using Git.

```bash
git clone git@github.com:Cinnamonness/BioFilterToolsPro.git
cd BioFilterToolsPro
```

---

## Features

### Features of BioFilterToolsPro

The utility provides two main functions: 

1.  **Operations on DNA and RNA sequences:**
    - Determination of the type of molecule (DNA or RNA)
    - Transcription
    - Reverse order ("reversal")
    - Definition of a complementary sequence
    - Definition of a reverse complementary sequence

2. **Sequence filtering from FASTQ file based on:**
    - GC-content
    - Sequence length
    - Threshold value of the average reed quality (phred33 scale)

### Features of bio_files_processor

The utility provides three main functions: 

1. **Working with FASTA files:**
    - Converting a multi-line FASTA file to a single-line format
    - Obtaining a valid FASTA file

2. **Working with BLAST output files:**
    - Writing proteins with the best matches from the BLAST output file
    - Generating a .txt file with a list of proteins sorted in alphabetical order

3. **Working with genomic annotations in .gbk format:**
    - Writing a specified number of genes before and after the given genes of interest along with their protein sequences (translation)
    - Obtaining a valid FASTA file

---

## Requirements
- Python 3.x
- Required libraries you can see in requirements.txt

---

## Example: 

### Working example with DNA sequence:

```Python
dna = DNASequence("ATGCGA")
print(dna)
print(dna[1:4])
print(dna[1::])
print(dna.complement())
print(dna.reverse())
print(dna.reverse_complement())
print(dna.transcribe())
```
```Python
Sequence: ATGCGA
Sequence: TGC
Sequence: TGCGA
Sequence: TACGCT
Sequence: AGCGTA
Sequence: TCGCAT
Sequence: AUGCGA
```

### Working example with RNA sequence:
```Python
rna = RNASequence("AUGCGA")
print(str(rna))        
print(rna[1:4])
print(rna[1::])
print(rna.complement())
print(rna.reverse())
print(rna.reverse_complement())
```
```Python
Sequence: AUGCGA
Sequence: UGC
Sequence: UGCGA
Sequence: UACGCU
Sequence: AGCGUA
Sequence: UCGCAU
```

### Working example with protein sequence:
```Python
aa_seq = AminoAcidSequence("GALNQRHKTTYCC")
print(aa_seq)
print(aa_seq[1:4])
print(aa_seq[1::])
print(aa_seq.classify_aminoacids())
```
```Python
Sequence: GALNQRHKTTYCC
Sequence: ALN
Sequence: ALNQRHKTTYCC
non-polar: Count = 3, Percentage = 23.08%
polar uncharged: Count = 7, Percentage = 53.85%
polar negatively charged: Count = 0, Percentage = 0.00%
polar positively charged: Count = 3, Percentage = 23.08%
```

### Working example of `filter_fastq`
```Python
if __name__ == "__main__":
    arguments = ("../data/example.fastq", 
             output_fastq='../data/filtered.fastq',
             gc_bounds=(0, 80),
             length_bounds=(0, 500),
             quality_threshold=40)
    result = filter_fastq(*arguments)
```
As a result of running the filter_fastq function, a file named filtered.fastq containing the filtered sequences from the original file example.fastq will be saved in the ./data directory. If the filtered directory did not previously exist, it will be created in the current directory.

***Original example_fastq.fastq***
```Python
```Python
@K00271:89:HHWWNBBXX:2:1101:23277:1068 1:N:0:CAGATC
NATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAATCTCAGACAACAAATCACAGAGTTAAGTCAGTTTACCGCACAAACTNACCAAGTCGGCGAAACAGAAGGTGGCGAC
+
#AAAFFJJFJJFJJJJF-FFFJJFJJJJFJJJJFJJFAFJF-FFJFJFJJJJ7A<FFJJJFAJF<<-JJJFJJF----FA--7-A-7------7A-7----7FJ----7---7-7A-AF-#A<---7---7A-F)-7AA<--77-)--)-
@K00271:89:HHWWNBBXX:2:1101:16792:1121 1:N:0:CAGATC
NGGAAACCGTCGGTTCTGGTTGTGGAGGCGGTTGGTGGTGGCTGTGGATTGGGAAGATGGTTTGGTAAGTTTTGGGCCGGTTTCAAGAAACTAAGCTGGGCTGGGCTGGGCTGGGCTGGGCTAAGCTGGGCTCACGAACCAAGTAAAGTT
+
#AAFFJJJJJJJJJJJJJJFFJFJJJJJFJJ-FFJFFJFJJJJJJJJFFFJJJ<<JF<FJ<FJJJFFJJJJFJJJJJJJJ-AJJFJJJAJJFJJJJFJJJJJJFFJJJJJJFJJJJJJJJJ-AFJFJJJJJFFAJJ-AFJJJFF-7F<77
@K00271:89:HHWWNBBXX:2:1101:18609:1173 1:N:0:CAGATC
NTGAACCTCCAAATTGATTAGAGAGTCACATGCAAATCACTGTTAAAGTCCGACAAACAGTTTACAAGTAACTTTAATTGAATGCCCGAAAATCTTTTAAATGTTCAGACAATACGATGACGATGACAGTTATAAGCGAAACCACTGAGA
+
#AAFFJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ<FJJJFJFJJJJJJJFJJJJFJJFJJJF
@K00271:89:HHWWNBBXX:2:1101:19898:1226 1:N:0:CAGATC
NTGATGACTGTTGCCAAACAATTTGGGAATTCTAGATGGGATTCGAGTTTAGTTTTGGAGTGAGCCTTATAATTTTGGTTCATCAAGGTCAATAAGGATACACTCCCACATTGGTGTTCATTGGGTTAATTTTGGAGTGCCACTCACACC
+
#AAFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ<FFJJFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJAJFJJJJJJJJJJJFJJFJJJJJJJJJJJJJJJJJJ
```
***Filtered filtered_sequences.fastq***
```Python
@K00271:89:HHWWNBBXX:2:1101:18609:1173 1:N:0:CAGATC
NTGAACCTCCAAATTGATTAGAGAGTCACATGCAAATCACTGTTAAAGTCCGACAAACAGTTTACAAGTAACTTTAATTGAATGCCCGAAAATCTTTTAAATGTTCAGACAATACGATGACGATGACAGTTATAAGCGAAACCACTGAGA
+
#AAFFJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ<FJJJFJFJJJJJJJFJJJJFJJFJJJF
@K00271:89:HHWWNBBXX:2:1101:19898:1226 1:N:0:CAGATC
NTGATGACTGTTGCCAAACAATTTGGGAATTCTAGATGGGATTCGAGTTTAGTTTTGGAGTGAGCCTTATAATTTTGGTTCATCAAGGTCAATAAGGATACACTCCCACATTGGTGTTCATTGGGTTAATTTTGGAGTGCCACTCACACC
+
#AAFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ<FFJJFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJAJFJJJJJJJJJJJFJJFJJJJJJJJJJJJJJJJJJ
```

### Working example of `convert_multiline_fasta_to_oneline` from ***bio_files_processor***

```Python
if __name__ == "__main__":
    arguments = ('C:/Users/User/Downloads/example_multiline_fasta.fasta', '')
    result = convert_multiline_fasta_to_oneline(*arguments)
```
***Original multi-line example_multiline_fasta.fasta***
```Python
>5S_rRNA::NODE_272_length_223_cov_0.720238:18-129(+)
ACGGCCATAGGACTTTGAAAGCACCGCATCCCGTCCGATCTGCGAAGTTAACCAAGATGCCGCCTGGTTAGTACCATGGTGGGGGACCACATGGGAATCCCT
GGTGCTGTG
>16S_rRNA::NODE_4_length_428221_cov_75.638017:281055-282593(-)
TTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAACAGCTTGCTGTTTCGCTGACGAGTGG
CGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACC
TTCGGGCCTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTTGTTGGTGGGGTAACGGCTCACCAAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCA
GCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGA
AGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTC
CGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCG
GGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAAT
ACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGA
TGTCGACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATT
GACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACGGAAGTTTTCAGAGATGAGAAT
GTGCCTTCGGGAACCGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGT
TGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACA
CACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCA
TGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAA
GAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGGGAACCTGCGGTTGGATC
ACCTCCTT
```
***One-line output_fasta.fasta***
```Python
>5S_rRNA::NODE_272_length_223_cov_0.720238:18-129(+)
ACGGCCATAGGACTTTGAAAGCACCGCATCCCGTCCGATCTGCGAAGTTAACCAAGATGCCGCCTGGTTAGTACCATGGTGGGGGACCACATGGGAATCCCTGGTGCTGTG
>16S_rRNA::NODE_4_length_428221_cov_75.638017:281055-282593(-)
TTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAACAGCTTGCTGTTTCGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTTGTTGGTGGGGTAACGGCTCACCAAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACGGAAGTTTTCAGAGATGAGAATGTGCCTTCGGGAACCGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGGGAACCTGCGGTTGGATCACCTCCTT
```

### Working example of `parse_blast_output` from ***bio_files_processor***

```Python
if __name__ == "__main__":
    arguments = ('C:/Users/User/Downloads/example_blast_results.txt', '')
    result = parse_blast_output(*arguments)
```

You can find the results of the parse_blast_output function in the examples directory. example_blast_results.txt is the original file that needed to be parsed. parse_blast.txt is the file containing the selected proteins.

### Working example of `select_genes_from_gbk_to_fasta` from ***bio_files_processor***
```Python
if __name__ == "__main__":
    arguments = ('C:/Users/User/Downloads/example_gbk.gbk', 'sucA', 'btuD_1')
    result = select_genes_from_gbk_to_fasta(*arguments, n_before=2, n_after=2, output_fasta='')
```

The results of the select_genes_from_gbk_to_fasta function can be found in the examples directory. example_gbk.gbk is the original file containing the genome annotation of ***E. coli*** gbk.fasta is the FASTA file with the selected number of genes before and after each gene of interest, along with the corresponding protein sequences (translation).
