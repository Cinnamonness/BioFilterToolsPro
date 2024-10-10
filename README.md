# BioFilterToolsPro

-[russian version](#russian-version)

-[english version](#english-version)

### russian version

**BioFilterToolsPro** — это утилита, предназначенная для работы с последовательностями ДНК и РНК, а также для фильтрации последовательностей FASTQ-файла на основе GC-состава, длины рида и порогового значения среднего качества рида (шкала phred33). 

Авторы:
* **Программное обеспечение:** *Карицкая Полина* <cinnamonness@gmail.com>, <br/>
Институт Биоинформатики, Санкт-Петербург, Россия. 
* **Идея, руководитель:** *Никита Ваулин*, *Антон Сидорин* <br/>
Институт Биоинформатики, Санкт-Петербург, Россия.

---
## Оглавление
- [Установка](#установка)
- [Особенности](#особенности)
- [Системные требования](#системные-требования)
- [Пример использования](#пример-использования-)

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
---
# BioFilterToolsPro

### english version

**BioFilterToolsPro** — is a utility designed to work with DNA and RNA sequences, as well as to filter sequences in FASTQ file based on GC composition, reed length and the threshold value of the average reed quality (phred33 scale).

Authors:
* **Software:** *Karitskaya Polina* <cinnamonness@gmail.com>, <br/>
Institute of Bioinformatics, Saint-Petersburg, Russia. 
* **Idea, supervisor:** *Nikita Vaulin*, *Anton Sidorin* <br/>
Institute of Bioinformatics, Saint-Petersburg, Russia.

---
## Table of contents
- [How to install](#how-to-install)
- [Features](#features)
- [Requirements](#requirements)
- [Example](#example-)

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
