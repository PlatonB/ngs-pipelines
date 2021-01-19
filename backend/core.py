__version__ = 'v3.0'

import os, re, sys
from subprocess import run, PIPE, STDOUT

class DifFmtsError(Exception):
        '''
        В исходной папке должны быть только одноформатные файлы.
        '''
        def __init__(self, src_file_fmts):
                err_msg = f'\nИсходные файлы разных форматов: {src_file_fmts}'
                super().__init__(err_msg)
                
class UnsupFmtError(Exception):
        '''
        Если папка загрязнена файлами с явно нефастовыми
        расширениями, скорее всего, исследователь не
        разобрался в предназначении проги. В таком
        случае, безопаснее будет выполнение прервать.
        '''
        def __init__(self, src_file_fmt):
                err_msg = f'\n{src_file_fmt} - неподдерживаемый формат исходного файла'
                super().__init__(err_msg)
                
class UnevenPairedReadsQuanError(Exception):
        '''
        То неловкое чувство,
        когда у парноридовой фасты
        нет второй половинки.
        '''
        def __init__(self):
                err_msg = '''\nПри исследовании парных ридов количество
исходных файлов должно быть чётным'''
                super().__init__(err_msg)
                
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
                Костяк атрибутов создаётся на основе общих
                аргументов. Фронтенды могут потом создать
                дополнительные атрибуты по своим аргументам.
                '''
                self.src_dir_path = os.path.normpath(args.src_dir_path)
                self.src_file_names = os.listdir(self.src_dir_path)
                src_file_fmts = set(map(lambda src_file_name:
                                        src_file_name.replace('.gz', '').rsplit('.', maxsplit=1)[1],
                                        self.src_file_names))
                if len(src_file_fmts) > 1:
                        raise DifFmtsError(src_file_fmts)
                self.src_file_fmt = list(src_file_fmts)[0]
                if self.src_file_fmt not in ['fq', 'fastq', 'fa', 'fasta', 'mfa', 'fna']:
                        raise UnsupFmtError(self.src_file_fmt)
                self.ref_file_path = os.path.normpath(args.ref_file_path)
                self.ref_dir_path, self.ref_file_name = os.path.split(self.ref_file_path)
                self.species_name = re.sub(r'[\.-]', '_', args.species_name).lower()
                self.reads_type = args.reads_type
                if args.trg_top_dir_path is None:
                        self.trg_top_dir_path = self.src_dir_path
                else:
                        self.trg_top_dir_path = os.path.normpath(args.trg_top_dir_path)
                self.threads_quan = args.threads_quan
                
        def compress_ref_file(self):
                '''
                Если референсный геном изначально не заархивирован
                с помощью BGZip, то это будет осуществлено при
                вызове данной функции. Но имейте в виду, что
                конвертация GZip в BGZip пока не поддерживается.
                '''
                if self.ref_file_name[-3:] != '.gz':
                        print(f'\n{self.ref_file_name}: сжатие\n')
                        run_command(f'bgzip -@ {self.threads_quan} -l 9 -i {self.ref_file_path}')
                        self.ref_file_name += '.gz'
                        self.ref_file_path = os.path.join(self.ref_dir_path,
                                                          self.ref_file_name)
                        
        def index_ref_file(self):
                '''
                Для многих задач (в частности, выравнивания и коллинга)
                референсный геном требуется проиндексировать с помощью
                bowtie2-build и samtools faidx. Если это было сделано
                уже раньше, то данные действия будут пропущены.
                '''
                for file_name in os.listdir(self.ref_dir_path):
                        if re.match(self.species_name + r'\.\d+\.bt2', file_name) is not None:
                                break
                else:
                        print(f'\n{self.ref_file_name}: создание .bt2-индексов\n')
                        run_command(f'bowtie2-build --threads {self.threads_quan} {self.ref_file_path} {os.path.join(self.ref_dir_path, self.species_name)}')
                if not os.path.exists(f'{self.ref_file_path}.fai'):
                        print(f'\n{self.ref_file_name}: создание .fai-индекса\n')
                        run_command(f'samtools faidx {self.ref_file_path}')
                        
        def group_file_names(self):
                '''
                Риды могут быть парными или непарными.
                Парные риды поставляются в виде пар
                файлов, что осложняет перебор содержимого
                исходной папки. Функция упорядочивает
                имена файлов, требуя лишь минимального
                вмешательства исследователя: ему нужно
                при запуске программы указывать тип ридов.
                '''
                
                #Сортировка имён исследуемых FASTA/Q. В
                #случае, если риды парные, то в получаемом
                #списке файлы каждого сэмпла оказываются
                #друг относительно друга по соседству.
                src_file_names = sorted(self.src_file_names)
                
                #Пайплайны предполагается применять отдельно для каждого
                #(или каждой пары) FASTA/Q. В качестве подготовки к
                #считыванию исходных данных растащим имена файлов,
                #будь они парные или непарные, по вложенным спискам.
                src_files_quan = len(src_file_names)
                if self.reads_type == 'paired':
                        if src_files_quan % 2 == 1:
                                raise UnevenPairedReadsQuanError()
                        else:
                                src_file_grps = [src_file_names[index:index + 2] for index in range(0, src_files_quan, 2)]
                else:
                        src_file_grps = [[src_file_name] for src_file_name in src_file_names]
                        
                return src_file_grps
        
        def get_srtd_bam(self, src_files_grp):
                '''
                Выравнивание ридов на референс
                и получение максимально готового
                к дальнейшей обработке BAM-файла.
                '''
                
                #Bowtie 2 требует точного определения формата
                #выравниваемых файлов: FASTA это или FASTQ.
                #Программа распознает формат по расширению.
                if self.src_file_fmt in ['fq', 'fastq']:
                        reads_fmt_opt = '-q'
                elif self.src_file_fmt in ['fa', 'fasta', 'mfa', 'fna']:
                        reads_fmt_opt = '-f'
                        
                #Базовая часть имени файла с результатами обработки пары
                #перекочует из общей для этой пары начальной части имён. Если
                #FASTA/Q непарные, то базовая часть имени каждого конечного файла
                #будет представлять собой имя соответствующего FASTA/Q без расширений.
                if self.reads_type == 'paired':
                        src_file_1_path = os.path.join(self.src_dir_path,
                                                       src_files_grp[0])
                        src_file_2_path = os.path.join(self.src_dir_path,
                                                       src_files_grp[1])
                        reads_files_opt = f'-1 {src_file_1_path} -2 {src_file_2_path}'
                        min_letters_quan = min(len(src_files_grp[0]),
                                               len(src_files_grp[1]))
                        src_file_base = ''
                        for index in range(min_letters_quan):
                                if src_files_grp[0][index] == src_files_grp[1][index]:
                                        src_file_base += src_files_grp[0][index]
                                else:
                                        break
                else:
                        src_file_path = os.path.join(self.src_dir_path,
                                                     src_files_grp[0])
                        reads_files_opt = f'-U {src_file_path}'
                        src_file_base = src_files_grp[0].replace('.gz', '').replace(f'.{self.src_file_fmt}', '')
                        
                #1 исследуемый FASTA/Q (или 1 пара) - 1 конечная подпапка.
                trg_dir_path = os.path.join(self.trg_top_dir_path,
                                            src_file_base)
                os.mkdir(trg_dir_path)
                
                #Создание пути к будущему BAM-файлу - продукту
                #выравнивания и трёх этапов SAMtools-обработки.
                #Путь будет в двух вариантах - без суффикса и
                #расширения, и с таковыми. Первый может пригодиться
                #для наименования продуктов пост-BAM-звеньев того
                #или иного пайплайна. Второй нужен в этой функции.
                bam_file_base_path = os.path.join(trg_dir_path,
                                                  src_file_base)
                bam_file_path = bam_file_base_path + '_srtd.bam'
                
                #Универсальная часть большого разнообразия возможных пайплайнов. Выравнивание
                #исследуемого FASTA/Q на референсный геном. SAM-вывод сразу перехватывается
                #для получения и сортировки BAM. Отсортированный BAM индексируется.
                print(f'\n{", ".join(src_files_grp)}: выравнивание, SAM --> BAM, сортировка и индексация BAM\n')
                run_command(f'''
bowtie2 -p {self.threads_quan} -x {os.path.join(self.ref_dir_path, self.species_name)} {reads_fmt_opt} {reads_files_opt} |
samtools view -@ {self.threads_quan} -b |
samtools sort -@ {self.threads_quan} -o {bam_file_path} &&
samtools index -@ {self.threads_quan} {bam_file_path}''')
                
                return bam_file_base_path
