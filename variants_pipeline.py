__version__ = 'v3.2'

import sys, os, re, gzip
sys.dont_write_bytecode = True
from cli.ngs_pipelines_cli_ru import add_args_ru
from backend.core import Core, run_command
from backend.prep_dbsnp_data import prep_dbsnp_data
from contextlib import ExitStack
from pysam import TabixFile, asTuple

tool_descr = '''
Python3-программа, формирующая пайплайн
от ридов до коротких герминальных вариантов.

Упрощённая схема пайплайна:
Выравнивание --> коллинг --> (только для человека) поиск rsIDs
'''

#Добавление кросспайплайновых аргументов.
args = add_args_ru(tool_descr, __version__)

#Создание объекта класса, содержащего
#общие для разнообразных пайплайнов
#методы. Добавление атрибута,
#обозначающего тип секвенирования.
core = Core(args)
core.seq_type = args.seq_type.upper()

#Архивация (BGZIP) и индексация
#(Bowtie2 и Samtools) FASTA/Q
#референсного генома, скачивание
#данных dbSNP и формирование их
#урезанной версии, группировка
#имён исследуемых FASTA/Q.
core.compress_ref_file()
core.index_ref_file()
if core.species_name == 'homo_sapiens':
        dbsnp_tsv_paths = prep_dbsnp_data(core.ref_dir_path,
                                          core.threads_quan)
src_file_grps = core.group_file_names()

#Перебор списков, в каждом из
#которых одно или два имени FASTA/Q.
for src_files_grp in src_file_grps:
        
        #Выравнивание, конвертация SAM в
        #BAM, сортировка и индексация BAM.
        bam_file_base_path = core.get_srtd_bam(src_files_grp)
        
        #Неинтересная возня с путями к файлам. Когда наладится
        #запуск DeepVariant без Docker, этот код упростится.
        trg_dir_path, bam_file_name = os.path.split(f'{bam_file_base_path}_srtd.bam')
        rawvcf_file_path = f'{bam_file_base_path}.vcf.gz'
        rawvcf_file_name = os.path.basename(rawvcf_file_path)
        vcf_file_path = f'{bam_file_base_path}_ann.vcf'
        
        #По BAM и референсному геному
        #осуществляется коллинг SNPs.
        print(f'\n{bam_file_name}: коллинг\n')
        run_command(f'''
docker run \
-v "{core.ref_dir_path}":"/ref" \
-v "{trg_dir_path}":"/trg" \
google/deepvariant /opt/deepvariant/bin/run_deepvariant \
--num_shards={core.threads_quan} --model_type={core.seq_type} \
--ref=/ref/{core.ref_file_name} \
--reads=/trg/{bam_file_name} \
--output_vcf=/trg/{rawvcf_file_name}''')
        
        #Ниже этого блока будет код,
        #применимый только к геному человека.
        if core.species_name != 'homo_sapiens':
                continue
        
        #Полученный в результате коллинга VCF не содержит
        #идентификаторов вариантов. В случае человеческого генома
        #их можно получить по dbSNP, что будет сделано ниже.
        print(f'\n{rawvcf_file_name}: фильтрация и замена точек на rsIDs\n')
        with gzip.open(rawvcf_file_path, mode='rt') as rawvcf_file_opened:
                
                #Предполагается, что варианты dbSNP распределены по отдельным
                #файлам (1 архив - 1 хромосома). По этой причине для разных
                #аннотируемых SNPs нужно доставать rsIDs из разных dbSNP-архивов.
                #Многочисленные открытия-закрытия архивов займут много времени,
                #поэтому откроем на чтение парсером pysam их всех сразу.
                with ExitStack() as stack:
                        dbsnp_tsvs_opened = {os.path.basename(dbsnp_tsv_path)[:-7]:
                                             stack.enter_context(TabixFile(dbsnp_tsv_path)) for dbsnp_tsv_path in dbsnp_tsv_paths}
                        
                        #Открытие конечного VCF на чтение и перегонка
                        #в него хэдеров свежеколленного VCF.
                        with open(vcf_file_path, 'w') as vcf_file_opened:
                                for line in rawvcf_file_opened:
                                        if line.startswith('##'):
                                                vcf_file_opened.write(line)
                                        else:
                                                vcf_file_opened.write(line)
                                                break
                                        
                                #В dbSNP-VCF может найтись несколько соответствий координате текущего
                                #запрашиваемого варианта. Как правило, это когда для данной позиции известны
                                #несколько разных вставок. Ради выявления вхождения наколленного варианта
                                #в dbSNP-референс сверяются не только координаты, но и соответствующие
                                #аллели. Если вхождение обнаружилось по координате и подтвердилось
                                #по аллелям, точка из сырого VCF заменяется на rsID из dbSNP.
                                for line in rawvcf_file_opened:
                                        vcf_row = line.split('\t')
                                        chrom, pos, qual = vcf_row[0], int(vcf_row[1]), float(vcf_row[5])
                                        if qual < 20:
                                                continue
                                        try:
                                                for dbsnp_tup in dbsnp_tsvs_opened[chrom].fetch(chrom,
                                                                                                pos - 1,
                                                                                                pos,
                                                                                                parser=asTuple()):
                                                        if vcf_row[3] == dbsnp_tup[3] and vcf_row[4] in dbsnp_tup[4].split(','):
                                                                vcf_row[2] = dbsnp_tup[2]
                                                                vcf_file_opened.write('\t'.join(vcf_row))
                                                                break
                                                else:
                                                        vcf_file_opened.write(line)
                                        except KeyError:
                                                vcf_file_opened.write(line)
                                                
        #BGZIP-архивация окончательного VCF и
        #удаление вместе с индексом сырого VCF.
        print(f'\n{os.path.basename(vcf_file_path)}: сжатие')
        run_command(f'bgzip -@ {core.threads_quan} -l 9 -i {vcf_file_path}')
        os.remove(rawvcf_file_path)
        os.remove(rawvcf_file_path + '.tbi')
