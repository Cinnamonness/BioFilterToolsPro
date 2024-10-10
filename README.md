# BioFilterToolsPro

-[russian version](#russian-version)

-[english version](#english-version)

### russian version

**BioFilterToolsPro** — это утилита, предназначенная для работы с последовательностями ДНК и РНК, а также для фильтрации последовательностей FASTQ-файла на основе GC-состава, длины рида и порогового значения среднего качества рида (шкала phred33). 

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
Необходимые модули для работы (`module_rna_dna_tools`, `module_filter_fastq`) находятся внутри директории этого проекта. 

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
- Необходимые модули:
    - `module_rna_dna_tools`
    - `module_filter_fastq`

---

## Пример использования: 

### Пример работы `run_dna_rna_tools`
```Python
if __name__ == "__main__":
    arguments = ('ATG', 'aT', 'reverse')
    result = run_dna_rna_tools(*arguments)
    print("Result of run_dna_rna_tools:", result)
```
```Python
Result of run_dna_rna_tools: ['GTA', 'Ta']
```
### Пример работы `filter_fastq`
```Python
if __name__ == "__main__":
    arguments = ('C:/Users/User/PycharmProjects/some_scripts/example_fastq.fastq', '',
                 (50, 70), (20,50), 20)
    result = filter_fastq(*arguments)
```
В результате работы функции `filter_fastq` в директории ./filtered сохранится файл filtered_sequences.fastq с отфильтрованными последовательностями из изначального файла example_fastq.fastq.
Если директории filtered ранее не существовала, то она будет создана в текущей директории. 

***Исходный example_fastq.fastq***
```Python
@SRX079804:1:SRR292678:1:1101:278698:278698 1:N:0:1 BH:ok
CTAATAATGGTAATTGAACCATAGAAGATAAGTTCATAATGTAATAAATACATCCATAGAGTTATTAA
+SRX079804:1:SRR292678:1:1101:278698:278698 1:N:0:1 BH:ok
DDBDBCCCDD@FFFB9<<<@DA=DA@B:@=@@AC@GGFCGECFFDGGCGFFGGFFCEBF9>?@>BDFF
@SRX079804:1:SRR292678:1:1101:918742:918742 1:N:0:1 BH:failed
CTCTCCATGCACAAAGAATATCACAGCCAAA
+SRX079804:1:SRR292678:1:1101:918742:918742 1:N:0:1 BH:failed
EEEBA?@;B@EEE@BEE=?EDDDDADCDA?E
@SRX079804:1:SRR292678:1:1101:923787:923787 2:N:0:1 BH:ok
TTGTGAAGGATGGGATATTAGTGTAGATGA
+SRX079804:1:SRR292678:1:1101:923787:923787 2:N:0:1 BH:ok
EEBBEGEEE=BBB<@DCDCGD@D>=DEGEE
@SRX079804:1:SRR292678:1:1101:933189:933189 1:N:0:1 BH:failed
GTCTGCACTATCGAGGGCTGTGCCTTTGC
+SRX079804:1:SRR292678:1:1101:933189:933189 1:N:0:1 BH:failed
FEFFDBFF8FE>?DFFFCEBCEEBBEDE6
@SRX079804:1:SRR292678:1:1101:937136:937136 1:N:0:1 BH:failed
TTTCTTTGGCTTAAAGATAGTTTTAGTC
+SRX079804:1:SRR292678:1:1101:937136:937136 1:N:0:1 BH:failed
EFFFEEEEFCBCDDDDE@/E?@@7@@3<
@SRX079804:1:SRR292678:1:1101:940351:940351 1:N:0:1 BH:changed:1
TGCCGTGGGAATGACAAACAAGCATCC
+SRX079804:1:SRR292678:1:1101:940351:940351 1:N:0:1 BH:changed:1
DECC@GFFBF=EBEAFDFGD?FFF8FF
@SRX079804:1:SRR292678:1:1101:940693:940693 1:N:0:1 BH:failed
CACATTATGAACTATGGGCACTGCAT
+SRX079804:1:SRR292678:1:1101:940693:940693 1:N:0:1 BH:failed
EEEGFDEDFEGGGGGFEGBGGGFGGG
@SRX079804:1:SRR292678:1:1101:955819:955819 1:N:0:1 BH:failed
CACCTAGCAGCAACGGACGAGTCAG
+SRX079804:1:SRR292678:1:1101:955819:955819 1:N:0:1 BH:failed
GGGGGEEEGGEGGGFGEGG;F@EFF
@SRX079804:1:SRR292678:1:1101:958051:958051 2:N:0:1 BH:ok
TTAATATTTCCATCTGAACTTCGC
+SRX079804:1:SRR292678:1:1101:958051:958051 2:N:0:1 BH:ok
EDDBGFEGFGHHFHGGEDEGBGDB
@SRX079804:1:SRR292678:1:1101:996098:996098 1:N:0:1 BH:failed
CTAAGAGAGTTTGTAATGCGGAC
+SRX079804:1:SRR292678:1:1101:996098:996098 1:N:0:1 BH:failed
DD=DBDBDC4EFFFD@?CD@ACD
@SRX079804:1:SRR292678:1:1101:1020278:1020278 2:N:0:1 BH:ok
AAAGTGCAGAACATGCAGATAT
+SRX079804:1:SRR292678:1:1101:1020278:1020278 2:N:0:1 BH:ok
D>AC?GDDCD?DDADE@GABDG
```
***Отфильтрованный filtered_sequences.fastq***
```Python
@SRX079804:1:SRR292678:1:1101:933189:933189 1:N:0:1 BH:failed
GTCTGCACTATCGAGGGCTGTGCCTTTGC
+SRX079804:1:SRR292678:1:1101:933189:933189 1:N:0:1 BH:failed
FEFFDBFF8FE>?DFFFCEBCEEBBEDE6
@SRX079804:1:SRR292678:1:1101:940351:940351 1:N:0:1 BH:changed:1
TGCCGTGGGAATGACAAACAAGCATCC
+SRX079804:1:SRR292678:1:1101:940351:940351 1:N:0:1 BH:changed:1
DECC@GFFBF=EBEAFDFGD?FFF8FF
@SRX079804:1:SRR292678:1:1101:955819:955819 1:N:0:1 BH:failed
CACCTAGCAGCAACGGACGAGTCAG
+SRX079804:1:SRR292678:1:1101:955819:955819 1:N:0:1 BH:failed
GGGGGEEEGGEGGGFGEGG;F@EFF
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

**BioFilterToolsPro** — is a utility designed to work with DNA and RNA sequences, as well as to filter sequences in FASTQ file based on GC composition, reed length and the threshold value of the average reed quality (phred33 scale).

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
The necessary modules for operation (`module_rna_dna_tools`, `module_filter_fastq`) are located inside the directory of this project.

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
- Required modules:
    - `module_rna_dna_tools`
    - `module_filter_fastq`

---

## Example: 

### Working example of `run_dna_rna_tools`

if __name__ == "__main__":
    arguments = ('ATG', 'aT', 'reverse')
    result = run_dna_rna_tools(*arguments)
    print("Result of run_dna_rna_tools:", result)
```
```Python
Result of run_dna_rna_tools: ['GTA', 'Ta']
```

### Working example of `filter_fastq`

```Python
if __name__ == "__main__":
    arguments = ('C:/Users/User/PycharmProjects/some_scripts/example_fastq.fastq', '',
                 (50, 70), (20,50), 20)
    result = filter_fastq(*arguments)
```
As a result of running the filter_fastq function, a file named filtered_sequences.fastq containing the filtered sequences from the original file example_fastq.fastq will be saved in the ./filtered directory. If the filtered directory did not previously exist, it will be created in the current directory.

***Original example_fastq.fastq***
```Python
@SRX079804:1:SRR292678:1:1101:278698:278698 1:N:0:1 BH:ok
CTAATAATGGTAATTGAACCATAGAAGATAAGTTCATAATGTAATAAATACATCCATAGAGTTATTAA
+SRX079804:1:SRR292678:1:1101:278698:278698 1:N:0:1 BH:ok
DDBDBCCCDD@FFFB9<<<@DA=DA@B:@=@@AC@GGFCGECFFDGGCGFFGGFFCEBF9>?@>BDFF
@SRX079804:1:SRR292678:1:1101:918742:918742 1:N:0:1 BH:failed
CTCTCCATGCACAAAGAATATCACAGCCAAA
+SRX079804:1:SRR292678:1:1101:918742:918742 1:N:0:1 BH:failed
EEEBA?@;B@EEE@BEE=?EDDDDADCDA?E
@SRX079804:1:SRR292678:1:1101:923787:923787 2:N:0:1 BH:ok
TTGTGAAGGATGGGATATTAGTGTAGATGA
+SRX079804:1:SRR292678:1:1101:923787:923787 2:N:0:1 BH:ok
EEBBEGEEE=BBB<@DCDCGD@D>=DEGEE
@SRX079804:1:SRR292678:1:1101:933189:933189 1:N:0:1 BH:failed
GTCTGCACTATCGAGGGCTGTGCCTTTGC
+SRX079804:1:SRR292678:1:1101:933189:933189 1:N:0:1 BH:failed
FEFFDBFF8FE>?DFFFCEBCEEBBEDE6
@SRX079804:1:SRR292678:1:1101:937136:937136 1:N:0:1 BH:failed
TTTCTTTGGCTTAAAGATAGTTTTAGTC
+SRX079804:1:SRR292678:1:1101:937136:937136 1:N:0:1 BH:failed
EFFFEEEEFCBCDDDDE@/E?@@7@@3<
@SRX079804:1:SRR292678:1:1101:940351:940351 1:N:0:1 BH:changed:1
TGCCGTGGGAATGACAAACAAGCATCC
+SRX079804:1:SRR292678:1:1101:940351:940351 1:N:0:1 BH:changed:1
DECC@GFFBF=EBEAFDFGD?FFF8FF
@SRX079804:1:SRR292678:1:1101:940693:940693 1:N:0:1 BH:failed
CACATTATGAACTATGGGCACTGCAT
+SRX079804:1:SRR292678:1:1101:940693:940693 1:N:0:1 BH:failed
EEEGFDEDFEGGGGGFEGBGGGFGGG
@SRX079804:1:SRR292678:1:1101:955819:955819 1:N:0:1 BH:failed
CACCTAGCAGCAACGGACGAGTCAG
+SRX079804:1:SRR292678:1:1101:955819:955819 1:N:0:1 BH:failed
GGGGGEEEGGEGGGFGEGG;F@EFF
@SRX079804:1:SRR292678:1:1101:958051:958051 2:N:0:1 BH:ok
TTAATATTTCCATCTGAACTTCGC
+SRX079804:1:SRR292678:1:1101:958051:958051 2:N:0:1 BH:ok
EDDBGFEGFGHHFHGGEDEGBGDB
@SRX079804:1:SRR292678:1:1101:996098:996098 1:N:0:1 BH:failed
CTAAGAGAGTTTGTAATGCGGAC
+SRX079804:1:SRR292678:1:1101:996098:996098 1:N:0:1 BH:failed
DD=DBDBDC4EFFFD@?CD@ACD
@SRX079804:1:SRR292678:1:1101:1020278:1020278 2:N:0:1 BH:ok
AAAGTGCAGAACATGCAGATAT
+SRX079804:1:SRR292678:1:1101:1020278:1020278 2:N:0:1 BH:ok
D>AC?GDDCD?DDADE@GABDG
```
***Filtered filtered_sequences.fastq***
```Python
@SRX079804:1:SRR292678:1:1101:933189:933189 1:N:0:1 BH:failed
GTCTGCACTATCGAGGGCTGTGCCTTTGC
+SRX079804:1:SRR292678:1:1101:933189:933189 1:N:0:1 BH:failed
FEFFDBFF8FE>?DFFFCEBCEEBBEDE6
@SRX079804:1:SRR292678:1:1101:940351:940351 1:N:0:1 BH:changed:1
TGCCGTGGGAATGACAAACAAGCATCC
+SRX079804:1:SRR292678:1:1101:940351:940351 1:N:0:1 BH:changed:1
DECC@GFFBF=EBEAFDFGD?FFF8FF
@SRX079804:1:SRR292678:1:1101:955819:955819 1:N:0:1 BH:failed
CACCTAGCAGCAACGGACGAGTCAG
+SRX079804:1:SRR292678:1:1101:955819:955819 1:N:0:1 BH:failed
GGGGGEEEGGEGGGFGEGG;F@EFF
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
