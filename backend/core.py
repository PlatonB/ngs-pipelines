__version__ = 'V2.0'

import os, re, sys
from argparse import ArgumentParser
from subprocess import run, PIPE, STDOUT

def add_main_args():
        '''
        Эта функция, впрочем, как и
        весь модуль, основывается на
        той идее, что подавляющее число
        геномных пайплайнов имеет
        одну и ту же базовую часть.
        Последняя, по логике, должна
        предоставлять исследователю
        фиксированный набор параметров.
        Функция создаёт объект
        соответствующего парсера
        и наполняет его рядом
        основополагающих опций.
        Объект возвращается фронтендам,
        а они дальше могут добавить
        туда специфические опции.
        '''
        argparser = ArgumentParser(description='''
Краткая форма с большой буквы - обязательный аргумент.
В квадратных скобках - значение по умолчанию.
В фигурных скобках - перечисление возможных значений.
''')
        argparser.add_argument('-S', '--src-dir-path', metavar='str', dest='src_dir_path', type=str,
                               help='Путь к папке с исследуемыми FASTA/Q-файлами')
        argparser.add_argument('-G', '--ref-file-path', metavar='str', dest='ref_file_path', type=str,
                               help='Путь к FASTA/Q-файлу с референсным геномом')
        argparser.add_argument('-O', '--species-name', metavar='str', dest='species_name', type=str,
                               help='Биологический вид (допустимы точки, дефисы, заглавные и строчные буквы)')
        argparser.add_argument('-R', '--reads-type', dest='reads_type', type=str, choices=['paired', 'unpaired'],
                               help='Парные или непарные чтения')
        argparser.add_argument('-t', '--trg-top-dir-path', metavar='[None]', dest='trg_top_dir_path', type=str,
                               help='Путь к папке для результатов (по умолчанию - путь к исходной папке)')
        argparser.add_argument('-p', '--threads-quan', metavar='[4]', dest='threads_quan', type=str, default='4',
                               help='Количество потоков, задействуемых компонентами пайплайна (не переборщите ☺)')
        return argparser

'''
Запуск любой линуксовой команды
с выводом STDOUT и STDERR.
'''
run_command = lambda command: print(run(command,
                                        stdout=PIPE,
                                        stderr=STDOUT,
                                        universal_newlines=True,
                                        shell=True).stdout)

class Core():
        '''
        Класс содержит наиболее востребованные
        для обработки геномных данных функции.
        Если их всех поочерёдно вызвать, то
        можно получить отсортированный BAM.
        '''
        
        def __init__(self, args):
                '''
                Костяк атрибутов создаётся
                на основе общих опций.
                Фронтенды могут потом
                создать дополнительные
                атрибуты по своим опциям.
                '''
                self.src_dir_path = os.path.normpath(args.src_dir_path)
                self.ref_file_path = os.path.normpath(args.ref_file_path)
                self.ref_dir_path, self.ref_file_name = os.path.split(self.ref_file_path)
                self.species_name = re.sub(r'[\.-]', '_', args.species_name).lower()
                self.reads_type = args.reads_type
                if args.trg_top_dir_path == None:
                        self.trg_top_dir_path = self.src_dir_path
                else:
                        self.trg_top_dir_path = os.path.normpath(args.trg_top_dir_path)
                self.threads_quan = args.threads_quan
                
        def compress_ref_file(self):
                '''
                Если референсный геном изначально
                не заархивирован с помощью
                BGZip, то это будет осуществлено
                при вызове данной функции.
                Но имейте в виду, что
                конвертация GZip в BGZip
                пока не поддерживается.
                '''
                if self.ref_file_name[-3:] != '.gz':
                        run_command(f'bgzip -@ {self.threads_quan} -l 9 -i {self.ref_file_path}')
                        self.ref_file_name += '.gz'
                        self.ref_file_path = os.path.join(self.ref_dir_path, self.ref_file_name)
                        
        def index_ref_file(self):
                '''
                Для многих задач (в частности,
                выравнивания и коллинга) референсный
                геном требуется проиндекцировать с
                помощью bowtie2-build и samtools faidx.
                Если это было сделано уже раньше,
                то данные действия будут пропущены.
                '''
                for file_name in os.listdir(self.ref_dir_path):
                        if re.match(self.species_name + r'\.\d+\.bt\d+', file_name) != None:
                                break
                else:
                        run_command(f'bowtie2-build --threads {self.threads_quan} {self.ref_file_path} {os.path.join(self.ref_dir_path, self.species_name)}')
                if os.path.exists(f'{self.ref_file_path}.fai') == False:
                        run_command(f'samtools faidx {self.ref_file_path}')
                        
        def group_file_names(self):
                '''
                Риды могут быть парными или непарными.
                Парные риды поставляются в виде
                пар файлов, что осложняет перебор
                содержимого исходной папки.
                Функция упорядочивает имена файлов,
                требуя лишь минимального вмешательства
                исследователя: ему нужно при запуске
                программы указывать тип ридов.
                '''
                
                #Сортировка имён исследуемых FASTA/Q.
                #В случае, если риды парные, то в
                #получаемом списке файлы каждого
                #запуска секвенатора оказываются
                #друг относительно друга по соседству.
                src_file_names = sorted(os.listdir(self.src_dir_path))
                
                #Пайплайны предполагается применять отдельно
                #для каждого (или каждой пары) FASTA/Q.
                #В качестве подготовки к считыванию исходных
                #данных растащим имена файлов, будь они
                #парные или непарные, по вложенным спискам.
                src_file_names_quan = len(src_file_names)
                if self.reads_type == 'paired':
                        if src_file_names_quan % 2 == 1:
                                print(f'''При исследовании парных ридов количество
исходных файлов должно быть чётным''')
                                sys.exit()
                        else:
                                nested_file_names = [src_file_names[index:index + 2] for index in range(0, src_file_names_quan, 2)]
                else:
                        nested_file_names = [[src_file_name] for src_file_name in src_file_names]
                        
                return nested_file_names
        
        def get_srtd_bam(self, element):
                '''
                Выравнивание ридов на референс
                и получение максимально готового
                к дальнейшей обработке BAM-файла.
                '''
                
                #Bowtie 2 требует точного
                #определения формата выравниваемых
                #файлов: FASTA это или FASTQ.
                #Программа попытается распознать
                #формат по расширению.
                #В случае неуспеха файл или
                #пара будут проигнорированы.
                if re.search(r'\.fq\b|\.fastq\b', element[0]) != None:
                        reads_format_opt = '-q'
                elif re.search(r'\.fa\b|\.fasta\b|.mfa\b|.fna\b', element[0]) != None:
                        reads_format_opt = '-f'
                else:
                        print(f'{element[0].split(".")[-1]} - неподдерживаемое расширение')
                        return None
                
                #Базовая часть имени файла с результатами
                #обработки пары перекочует из общей
                #для этой пары начальной части имён.
                #Если FASTA/Q непарные, то базовая
                #часть имени каждого конечного файла
                #будет представлять собой имя
                #соответствующего FASTA/Q без расширения.
                if self.reads_type == 'paired':
                        src_file_1_path = os.path.join(self.src_dir_path, element[0])
                        src_file_2_path = os.path.join(self.src_dir_path, element[1])
                        reads_files_opt = f'-1 {src_file_1_path} -2 {src_file_2_path}'
                        min_letters_quan = min(len(element[0]), len(element[1]))
                        src_file_base = ''
                        for index in range(min_letters_quan):
                                if element[0][index] == element[1][index]:
                                        src_file_base += element[0][index]
                                else:
                                        break
                else:
                        src_file_path = os.path.join(self.src_dir_path, element[0])
                        reads_files_opt = f'-U {src_file_path}'
                        src_file_base = '.'.join(element[0].split('.')[:-1])
                        
                #1 исследуемый FASTA/Q (или 1
                #пара) - 1 конечная подпапка.
                trg_dir_path = os.path.join(self.trg_top_dir_path, src_file_base)
                os.mkdir(trg_dir_path)
                
                #Создание пути к будущему
                #BAM-файлу - продукту
                #выравнивания и трёх
                #этапов SAMtools-обработки.
                #Путь будет в двух вариантах -
                #без расширения и с таковым.
                #Первый может пригодиться
                #для наименования продуктов
                #пост-BAM-звеньев того
                #или иного пайплайна.
                bam_file_base_path = os.path.join(trg_dir_path, src_file_base)
                bam_file_path = bam_file_base_path + '_srtd.bam'
                
                #Универсальная часть большого
                #разнообразия возможных пайплайнов.
                #Выравнивание исследуемого
                #FASTA/Q на референсный геном.
                #SAM-вывод сразу перехватывается
                #для получения и сортировки BAM.
                #Отсортированный BAM индексируется.
                run_command(f'''
bowtie2 -p {self.threads_quan} -x {os.path.join(self.ref_dir_path, self.species_name)} {reads_format_opt} {reads_files_opt} |
samtools view -@ {self.threads_quan} -b |
samtools sort -@ {self.threads_quan} -o {bam_file_path} &&
samtools index -@ {self.threads_quan} {bam_file_path}''')
                
                return bam_file_base_path
