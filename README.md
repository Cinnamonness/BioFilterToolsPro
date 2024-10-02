# BioFilterToolsPro

-[russian version](#russian-version)

-[english version](#english-version)

### russian version

**BioFilterToolsPro** — это утилита, предназначенная для работы с последовательностями ДНК и РНК, а также для фильтрации FASTQ-последовательностей на основе GC-состава, длины рида и порогового значения среднего качества рида (шкала phred33). 

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

2. **Фильтрация FASTQ-последовательности на основе:**
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

```Python
if __name__ == "__main__":
    arguments = ('ATG', 'aT', 'reverse')
    result = run_dna_rna_tools(*arguments)
    print("Result of run_dna_rna_tools:", result)
    seqs = {'@SRX079801':
            ('ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA',
             'FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD'),
    '@SRX079802':
        ('ATTAGCGAGGAGGAGTGCTGAGAAGATGTCGCCTACGCCGTTGAAATTCCCTTCAATCAGGGGGTACTGGAGGATACGAGTTTGTGTG',
         'BFFFFFFFB@B@A<@D>BDDACDDDEBEDEFFFBFFFEFFDFFF=CC@DDFD8FFFFFFF8/+.2,@7<<:?B/:<><-><@.A*C>D'),
    '@SRX079803':
                ('GAACGACAGCAGCTCCTGCATAACCGCGTCCTTCTTCTTTAGCGTTGTGCAAAGCATGTTTTGTATTACGGGCATCTCGAGCGAATC',
                 'DFFFEGDGGGGFGGEDCCDCEFFFFCCCCCB>CEBFGFBGGG?DE=:6@=>A<A>D?D8DCEE:>EEABE5D@5:DDCA;EEE-DCD'),
    '@SRX079804':
                ('TGAAGCGTCGATAGAAGTTAGCAAACCCGCGGAACTTCCGTACATCAGACACATTCCGGGGGGTGGGCCAATCCATGATGCCTTTG',
                 'FF@FFBEEEEFFEFFD@EDEFFB=DFEEFFFE8FFE8EEDBFDFEEBE+E<C<C@FFFFF;;338<??D:@=DD:8DDDD@EE?EB'),
    '@SRX079805':
                ('TAGGGTTGTATTTGCAGATCCATGGCATGCCAAAAAGAACATCGTCCCGTCCAATATCTGCAACATACCAGTTGGTTGGTA',
                 '@;CBA=:@;@DBDCDEEE/EEEEEEF@>FBEEB=EFA>EEBD=DAEEEEB9)99>B99BC)@,@<9CDD=C,5;B::?@;A')
    }
    gc_bounds = (10, 70)
    length_bounds = (10, 90)
    quality_threshold = 20
    filtered_seqs = filter_fastq(seqs, gc_bounds, length_bounds, quality_threshold)
    print("Result of filter_fastq:", filtered_seqs)
```
```Python
Result of run_dna_rna_tools: ['GTA', 'Ta']
Result of filter_fastq: {'@SRX079801': ('ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA', 'FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD'), 
'@SRX079802': ('ATTAGCGAGGAGGAGTGCTGAGAAGATGTCGCCTACGCCGTTGAAATTCCCTTCAATCAGGGGGTACTGGAGGATACGAGTTTGTGTG', 'BFFFFFFFB@B@A<@D>BDDACDDDEBEDEFFFBFFFEFFDFFF=CC@DDFD8FFFFFFF8/+.2,@7<<:?B/:<><-><@.A*C>D'), 
'@SRX079803': ('GAACGACAGCAGCTCCTGCATAACCGCGTCCTTCTTCTTTAGCGTTGTGCAAAGCATGTTTTGTATTACGGGCATCTCGAGCGAATC', 
'DFFFEGDGGGGFGGEDCCDCEFFFFCCCCCB>CEBFGFBGGG?DE=:6@=>A<A>D?D8DCEE:>EEABE5D@5:DDCA;EEE-DCD'), 
'@SRX079804': ('TGAAGCGTCGATAGAAGTTAGCAAACCCGCGGAACTTCCGTACATCAGACACATTCCGGGGGGTGGGCCAATCCATGATGCCTTTG', 
'FF@FFBEEEEFFEFFD@EDEFFB=DFEEFFFE8FFE8EEDBFDFEEBE+E<C<C@FFFFF;;338<??D:@=DD:8DDDD@EE?EB'), 
'@SRX079805': ('TAGGGTTGTATTTGCAGATCCATGGCATGCCAAAAAGAACATCGTCCCGTCCAATATCTGCAACATACCAGTTGGTTGGTA', 
'@;CBA=:@;@DBDCDEEE/EEEEEEF@>FBEEB=EFA>EEBD=DAEEEEB9)99>B99BC)@,@<9CDD=C,5;B::?@;A')}

```

---
# BioFilterToolsPro

### english version

**BioFilterToolsPro** — is a utility designed to work with DNA and RNA sequences, as well as to filter FASTQ sequences based on GC composition, reed length and the threshold value of the average reed quality (phred33 scale).

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

2. **FASTQ sequence filtering based on:**
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

```Python
if __name__ == "__main__":
    arguments = ('ATG', 'aT', 'reverse')
    result = run_dna_rna_tools(*arguments)
    print("Result of run_dna_rna_tools:", result)
    seqs = {'@SRX079801':
            ('ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA',
             'FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD'),
    '@SRX079802':
        ('ATTAGCGAGGAGGAGTGCTGAGAAGATGTCGCCTACGCCGTTGAAATTCCCTTCAATCAGGGGGTACTGGAGGATACGAGTTTGTGTG',
         'BFFFFFFFB@B@A<@D>BDDACDDDEBEDEFFFBFFFEFFDFFF=CC@DDFD8FFFFFFF8/+.2,@7<<:?B/:<><-><@.A*C>D'),
    '@SRX079803':
                ('GAACGACAGCAGCTCCTGCATAACCGCGTCCTTCTTCTTTAGCGTTGTGCAAAGCATGTTTTGTATTACGGGCATCTCGAGCGAATC',
                 'DFFFEGDGGGGFGGEDCCDCEFFFFCCCCCB>CEBFGFBGGG?DE=:6@=>A<A>D?D8DCEE:>EEABE5D@5:DDCA;EEE-DCD'),
    '@SRX079804':
                ('TGAAGCGTCGATAGAAGTTAGCAAACCCGCGGAACTTCCGTACATCAGACACATTCCGGGGGGTGGGCCAATCCATGATGCCTTTG',
                 'FF@FFBEEEEFFEFFD@EDEFFB=DFEEFFFE8FFE8EEDBFDFEEBE+E<C<C@FFFFF;;338<??D:@=DD:8DDDD@EE?EB'),
    '@SRX079805':
                ('TAGGGTTGTATTTGCAGATCCATGGCATGCCAAAAAGAACATCGTCCCGTCCAATATCTGCAACATACCAGTTGGTTGGTA',
                 '@;CBA=:@;@DBDCDEEE/EEEEEEF@>FBEEB=EFA>EEBD=DAEEEEB9)99>B99BC)@,@<9CDD=C,5;B::?@;A')
    }
    gc_bounds = (10, 70)
    length_bounds = (10, 90)
    quality_threshold = 20
    filtered_seqs = filter_fastq(seqs, gc_bounds, length_bounds, quality_threshold)
    print("Result of filter_fastq:", filtered_seqs)
```
```Python
Result of run_dna_rna_tools: ['GTA', 'Ta']
Result of filter_fastq: {'@SRX079801': ('ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA', 'FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD'), 
'@SRX079802': ('ATTAGCGAGGAGGAGTGCTGAGAAGATGTCGCCTACGCCGTTGAAATTCCCTTCAATCAGGGGGTACTGGAGGATACGAGTTTGTGTG', 'BFFFFFFFB@B@A<@D>BDDACDDDEBEDEFFFBFFFEFFDFFF=CC@DDFD8FFFFFFF8/+.2,@7<<:?B/:<><-><@.A*C>D'), 
'@SRX079803': ('GAACGACAGCAGCTCCTGCATAACCGCGTCCTTCTTCTTTAGCGTTGTGCAAAGCATGTTTTGTATTACGGGCATCTCGAGCGAATC', 
'DFFFEGDGGGGFGGEDCCDCEFFFFCCCCCB>CEBFGFBGGG?DE=:6@=>A<A>D?D8DCEE:>EEABE5D@5:DDCA;EEE-DCD'), 
'@SRX079804': ('TGAAGCGTCGATAGAAGTTAGCAAACCCGCGGAACTTCCGTACATCAGACACATTCCGGGGGGTGGGCCAATCCATGATGCCTTTG', 
'FF@FFBEEEEFFEFFD@EDEFFB=DFEEFFFE8FFE8EEDBFDFEEBE+E<C<C@FFFFF;;338<??D:@=DD:8DDDD@EE?EB'), 
'@SRX079805': ('TAGGGTTGTATTTGCAGATCCATGGCATGCCAAAAAGAACATCGTCCCGTCCAATATCTGCAACATACCAGTTGGTTGGTA', 
'@;CBA=:@;@DBDCDEEE/EEEEEEF@>FBEEB=EFA>EEBD=DAEEEEB9)99>B99BC)@,@<9CDD=C,5;B::?@;A')}

```
