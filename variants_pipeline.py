__version__ = 'V2.0'

print('''
Python3-программа, формирующая пайплайн
от ридов до характеристик SNPs.
Автор: Платон Быкадоров (platon.work@gmail.com), 2019-2020.
Версия: V2.0.
Лицензия: GNU General Public License version 3.
Поддержать проект: https://money.yandex.ru/to/41001832285976
Документация: https://github.com/PlatonB/ngs-pipelines/blob/master/README.md

Упрощённая схема пайплайна:
Выравнивание --> коллинг --> (только для человека) аннотирование.

На вход идут FASTA/Q-файлы с исследуемыми
ридами (WGS, WES, PacBio) и один
FASTA/Q с референсным геномом.
Каждый FASTA/Q-файл (или каждая пара
файлов) обрабатывается по-отдельности.

Если риды - парные, то каждый тандем исследуемых
файлов должен иметь общее начало имён:

SA1016N_R1.fastq.gz
SA1016N_R2.fastq.gz

SA1016T_R1.fastq.gz
SA1016T_R2.fastq.gz

SA1016X3_3L_R1.fastq.gz
SA1016X3_3L_R2.fastq.gz

Исследуемые риды и референсный геном могут
быть сжаты, но только с помощью BGZIP.
''')

import sys, os
sys.dont_write_bytecode = True
from backend.core import add_main_args, Core, run_command

#Добавление кросспайплайновых аргументов.
argparser = add_main_args()

#Специфическим аргументом
#этого пайплайна будет
#запрашиваемый программой
#DeepVariant тип секвенирования.
argparser.add_argument('--seq-type', dest='seq_type', type=str, choices=['wgs', 'wes', 'pacbio'],
                       help='Тип секвенирования')

#К парсингу разумно
#приступать только один
#раз, и притом после
#окончания формирования
#набора аргументов.
#Если парсить несколько
#раз, исследователю
#станет доступна только
#та часть аргументов,
#которая накопилась до
#момента первого парсинга.
args = argparser.parse_args()

#Создание объекта
#класса, содержащего
#общие для разнообразных
#пайплайнов методы.
#Добавление атрибута,
#обозначающего тип
#секвенирования.
core = Core(args)
core.seq_type = args.seq_type.upper()

#Архивация (BGZIP) и индексация
#(Bowtie2 и Samtools) FASTA/Q
#референсного генома, группирование
#имён исследуемых FASTA/Q.
core.compress_ref_file()
core.index_ref_file()
nested_file_names = core.group_file_names()

#Перебор списков, в каждом из
#которых одно или два имени FASTA/Q.
for element in nested_file_names:
        
        #Выравнивание, конвертация SAM в
        #BAM, сортировка и индексация BAM.
        bam_file_base_path = core.get_srtd_bam(element)
        
        #Неинтересная возня с путями к файлам.
        #Когда наладится запуск DeepVariant
        #без Docker, этот код упростится.
        if bam_file_base_path == None:
                continue
        trg_dir_path, bam_file_name = os.path.split(f'{bam_file_base_path}_srtd.bam')
        vcf_file_path = f'{bam_file_base_path}.vcf.gz'
        vcf_file_name = os.path.basename(vcf_file_path)
        ann_file_path = f'{bam_file_base_path}_ann.tsv'
        
        #По BAM и референсному геному
        #осуществляется коллинг SNPs.
        #SNPs человеческого генома
        #ещё и подробно аннотируются.
        run_command(f'''
sudo docker run \
-v "{core.ref_dir_path}":"/ref" \
-v "{trg_dir_path}":"/trg" \
google/deepvariant /opt/deepvariant/bin/run_deepvariant \
--num_shards={core.threads_quan} --model_type={core.seq_type} \
--ref=/ref/{core.ref_file_name} \
--reads=/trg/{bam_file_name} \
--output_vcf=/trg/{vcf_file_name}''')
        if core.species_name != 'homo_sapiens':
                continue
        run_command(f'''
vep --cache --fork {core.threads_quan} --quiet --no_stats --tab --offline -e \
-i {vcf_file_path} -o {ann_file_path}''')
        os.remove(vcf_file_path)
        os.remove(vcf_file_path + '.tbi')
