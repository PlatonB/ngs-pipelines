print('''
Python3-скрипт, формирующий пайплайн
от ридов до характеристик SNPs.
Автор: Платон Быкадоров (platon.work@gmail.com), 2019-2020.
Версия: V1.0.
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

def get_common_letters(file_name_1_char, file_name_2_char):
        '''
        Сопоставление букв имён исходных файлов.
        Нужно для выявления совпадающего начала.
        '''
        if file_name_1_char == file_name_2_char:
                return file_name_1_char
        
'''
Запуск любой линуксовой команды
с выводом STDOUT и STDERR.
'''
run_command = lambda command: print(run(command,
                                        stdout=PIPE,
                                        stderr=STDOUT,
                                        universal_newlines=True,
                                        shell=True).stdout)

####################################################################################################

import os, re, sys
from subprocess import run, PIPE, STDOUT

src_dir_path = os.path.normpath(input('Путь к папке с исследуемыми FASTA/Q-файлами: '))

ref_file_path = os.path.normpath(input('\nПуть к FASTA/Q-файлу с референсным геномом: '))
ref_dir_path, ref_file_name = os.path.split(ref_file_path)

trg_top_dir_path = os.path.normpath(input('\nПуть к папке для результатов: '))

species = re.sub(r'[\s\.-]', '_', input('\nБиологический вид: ')).lower()

reads = input('''\nЧтения - парные или непарные?
[paired(|p)|unpaired(|u)]: ''')
if reads not in ['paired', 'p', 'unpaired', 'u']:
        print(f'{reads} - недопустимая опция')
        sys.exit()
        
seq_type = input('''\nТип секвенирования
[wgs|wes|pacbio]: ''').upper()
if seq_type not in ['WGS', 'WES', 'PACBIO']:
        print(f'{seq_type} - недопустимая опция')
        sys.exit()
        
num_of_threads = input('''\nВо сколько потоков производить обработку?
(игнорирование ввода ==> 4 потока)
[1|2|3|4(|<enter>)|...]: ''')
if num_of_threads == '':
        num_of_threads = '4'
        
root_pass = input('''\nПароль root
(убедитесь, что никто не подглядывает ☺) ''')

#Если референсный геном изначально
#не заархивирован с помощью BGZip,
#то это будет осуществлено сейчас.
if ref_file_name[-3:] != '.gz':
        run_command(f'bgzip -@ 4 -l 9 -i {ref_file_path}')
        ref_file_name += '.gz'
        ref_file_path = os.path.join(ref_dir_path, ref_file_name)
        
#Референсный геном должен
#быть проиндексирован с
#помощью bowtie2-build
#(для дальнейшего выравнивания)
#и samtools faidx (для
#дальнейшего коллинга).
#Если это было сделано
#уже раньше, то данный
#этап будет пропущен.
for file_name in os.listdir(ref_dir_path):
        if re.match(species + r'\.\d+\.bt\d+', file_name) != None:
                break
else:
        run_command(f'bowtie2-build --threads {num_of_threads} {ref_file_path} {os.path.join(ref_dir_path, species)}')
if os.path.exists(f'{ref_file_path}.fai') == False:
        run_command(f'samtools faidx {ref_file_path}')
        
#Сортировка имён исследуемых FASTA/Q.
#В случае, если риды парные, то в
#получаемом списке файлы каждого
#запуска секвенатора оказываются
#друг относительно друга по соседству.
src_file_names = sorted(os.listdir(src_dir_path))

#Пайплайн должен будет применяться отдельно
#для каждого (или каждой пары) FASTA/Q.
#В качестве подготовки к считыванию исходных
#данных, растащим имена как парных, так
#и непарных файлов по вложенным спискам.
src_file_names_quan = len(src_file_names)
if reads in ['paired', 'p']:
        if src_file_names_quan % 2 == 1:
                print(f'''При исследовании парных ридов количество
исходных файлов должно быть чётным''')
                sys.exit()
        else:
                nested_file_names = [src_file_names[index:index + 2] for index in range(0, src_file_names_quan, 2)]
else:
        nested_file_names = [[src_file_name] for src_file_name in src_file_names]
        
#Перебор списков, в каждом из
#которых одно или два имени FASTA/Q.
for element in nested_file_names:
        
        #Bowtie 2 требует точного определения исходного
        #формата выравниваемого файла: FASTA это, или FASTQ.
        #Скрипт попытается распознать его по расширению.
        #В случае неуспеха файл будет проигнорирован.
        if re.search(r'\.fq\b|\.fastq\b', element[0]) != None:
                reads_format_opt = '-q'
        elif re.search(r'\.fa\b|\.fasta\b|.mfa\b|.fna\b', element[0]) != None:
                reads_format_opt = '-f'
        else:
                print(f'{src_file_name}: неподдерживаемое расширение ')
                continue
        
        #Базовая часть имени файла с результатами
        #обработки пары перекочует из общей
        #для этой пары начальной части имён.
        #Если FASTA/Q непарные, то базовая
        #часть имени каждого конечного файла
        #будет представлять собой имя
        #соответствующего FASTA/Q без расширения.
        if reads in ['paired', 'p']:
                src_file_1_path = os.path.join(src_dir_path, element[0])
                src_file_2_path = os.path.join(src_dir_path, element[1])
                reads_files_opt = f'-1 {src_file_1_path} -2 {src_file_2_path}'
                src_file_base_gen = map(get_common_letters, list(element[0]), list(element[1]))
                src_file_base = ''
                for char in src_file_base_gen:
                        if char != None:
                                src_file_base += char
                        else:
                                break
        else:
                src_file_path = os.path.join(src_dir_path, element[0])
                reads_files_opt = f'-U {src_file_path}'
                src_file_base = '.'.join(element[0].split('.')[:-1])
                
        #1 исследуемый FASTA/Q (или 1
        #пара) - 1 конечная подпапка.
        trg_dir_path = os.path.join(trg_top_dir_path, src_file_base)
        os.mkdir(trg_dir_path)
        
        #Создание путей к будущим результатам:
        #BAM-файлу (продукту выравнивания и
        #SAMtools-обработки), сжатому
        #VCF (продукту коллинга) и TSV
        #(получающемуся при аннотировании).
        bam_file_path = os.path.join(trg_dir_path, src_file_base) + '_srtd.bam'
        vcf_file_path = os.path.join(trg_dir_path, src_file_base) + '.vcf.gz'
        ann_file_path = os.path.join(trg_dir_path, src_file_base) + '_ann.tsv'
        
        #Основная часть пайплайна.
        #Выравнивание исследуемого
        #FASTA/Q на референсный геном.
        #SAM-вывод сразу перехватывается
        #для получения и сортировки BAM.
        #Отсортированный BAM индексируется.
        #По нему и референсному геному
        #осуществляется коллинг SNPs.
        #SNPs человеческого генома
        #ещё и подробно аннотируются.
        run_command(f'''
bowtie2 -p {num_of_threads} -x {os.path.join(ref_dir_path, species)} {reads_format_opt} {reads_files_opt} |
samtools view -@ {num_of_threads} -b |
samtools sort -@ {num_of_threads} -o {bam_file_path} &&
samtools index -@ {num_of_threads} {bam_file_path} &&
echo {root_pass} | sudo -S docker run \
-v "{ref_dir_path}":"/ref" \
-v "{trg_dir_path}":"/trg" \
google/deepvariant /opt/deepvariant/bin/run_deepvariant \
--num_shards={num_of_threads} --model_type={seq_type} \
--ref=/ref/{ref_file_name} \
--reads=/trg/{os.path.basename(bam_file_path)} \
--output_vcf=/trg/{src_file_base}.vcf.gz''')
        if species != 'homo_sapiens':
                continue
        run_command(f'''
vep --cache --fork {num_of_threads} --quiet --no_stats --tab --offline -e \
-i {vcf_file_path} -o {ann_file_path}''')
        os.remove(vcf_file_path)
        os.remove(vcf_file_path + '.tbi')
