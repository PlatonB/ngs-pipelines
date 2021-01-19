__version__ = 'v1.0'

from argparse import ArgumentParser, RawTextHelpFormatter

def add_args_ru(tool_descr, ver):
        '''
        Работа с аргументами командной строки.
        '''
        arg_parser = ArgumentParser(description=f'''
{tool_descr}

Версия: {ver}
Требуемые сторонние компоненты: Bowtie 2, SAMtools, DeepVariant
Автор: Платон Быкадоров (platon.work@gmail.com), 2019-2021
Лицензия: GNU General Public License version 3
Поддержать проект: https://www.tinkoff.ru/rm/bykadorov.platon1/7tX2Y99140/
Документация: https://github.com/PlatonB/ngs-pipelines/blob/master/README.md
Багрепорты/пожелания/общение: https://github.com/PlatonB/ngs-pipelines/issues

На вход идут FASTA/Q-файлы с исследуемыми
ридами (WGS, WES, PacBio) и один FASTA/Q с
референсным геномом. Каждый FASTA/Q-файл (или
каждая пара файлов) обрабатывается по-отдельности.

Исследуемые риды и референсный геном могут
быть сжаты, но только с помощью BGZIP.

Если риды - парные, то каждый тандем исследуемых
файлов должен иметь общее начало имён:

SA1016N_R1.fastq.gz
SA1016N_R2.fastq.gz

SA1016T_R1.fastq.gz
SA1016T_R2.fastq.gz

SA1016X3_3L_R1.fastq.gz
SA1016X3_3L_R2.fastq.gz

Условные обозначения в справке по CLI:
[значение по умолчанию на этапе парсинга аргументов];
[[конкретизированное значение по умолчанию]];
{{допустимые значения}}
''',
                                    formatter_class=RawTextHelpFormatter,
                                    add_help=False)
        hlp_grp = arg_parser.add_argument_group('Аргумент вывода справки')
        hlp_grp.add_argument('-h', '--help', action='help',
                             help='Вывести справку и выйти')
        man_grp = arg_parser.add_argument_group('Обязательные общие аргументы')
        man_grp.add_argument('-S', '--src-dir-path', required=True, metavar='str', dest='src_dir_path', type=str,
                             help='Путь к папке с исследуемыми FASTA/Q-файлами')
        man_grp.add_argument('-G', '--ref-file-path', required=True, metavar='str', dest='ref_file_path', type=str,
                             help='Путь к FASTA/Q-файлу с референсным геномом')
        man_grp.add_argument('-O', '--species-name', required=True, metavar='str', dest='species_name', type=str,
                             help='Биологический вид (допустимы точки, дефисы, заглавные и строчные буквы)')
        man_grp.add_argument('-R', '--reads-type', required=True, metavar='str', choices=['paired', 'unpaired'], dest='reads_type', type=str,
                             help='{paired, unpaired} Парные или непарные чтения')
        opt_grp = arg_parser.add_argument_group('Необязательные общие аргументы')
        opt_grp.add_argument('-t', '--trg-top-dir-path', metavar='[None]', dest='trg_top_dir_path', type=str,
                             help='Путь к папке для результатов ([[путь к исходной папке]])')
        opt_grp.add_argument('-p', '--threads-quan', metavar='[4]', default='4', dest='threads_quan', type=str,
                             help='Количество потоков, задействуемых компонентами пайплайна (не переборщите ☺)')
        gervar_grp = arg_parser.add_argument_group('Аргументы germline_variants_pipeline')
        gervar_grp.add_argument('-s', '--seq-type', metavar='[wgs]', choices=['wgs', 'wes', 'pacbio'], default='wgs', dest='seq_type', type=str,
                                help='{wgs, wes, pacbio} Тип секвенирования')
        args = arg_parser.parse_args()
        return args
