__version__ = 'V1.1'

import sys, os, urllib.request, re
sys.dont_write_bytecode = True
from backend.core import run_command

def prep_dbsnp_data(ref_dir_path, threads_quan):
        '''
        Функция соберёт из FTP проекта
        Ensembl таблицы со всеми известными
        генеративными мутациями (germline
        variants), срежет часть столбцов и
        проиндексирует полученный результат.
        Про срезание: целиком dbSNP-таблицы
        занимали бы многовато места, поэтому
        будут сохраняться только критически
        необходимые сниповому пайплайну
        первые 5 столбцов каждой из них.
        '''
        print('\ndbSNP: скачивание и обработка, если это не сделано ранее\n')
        
        #Сбор и размещение в файл
        #адресов архивов dbSNP.
        #Это действие, как
        #и все последующие в
        #функции, предполагается
        #выполнить однократно -
        #лишь при первом запуске
        #соответствующего пайплайна.
        dbsnp_urlstxt_path = os.path.join(ref_dir_path, 'dbsnp_urls.txt')
        ens_vcfspage_url = 'ftp://ftp.ensembl.org/pub/release-100/variation/vcf/homo_sapiens/'
        print('\ndbsnp_urls.txt', end='... ')
        if os.path.exists(dbsnp_urlstxt_path) == False:
                with urllib.request.urlopen(ens_vcfspage_url) as response:
                        dbsnp_vcf_names = re.findall(r'homo_sapiens-chr\w+\.vcf\.gz(?=\r\n)',
                                                     response.read().decode('UTF-8'))
                with open(dbsnp_urlstxt_path, 'w') as dbsnp_urlstxt_opened:
                        for dbsnp_vcf_name in dbsnp_vcf_names:
                                dbsnp_urlstxt_opened.write(os.path.join(ens_vcfspage_url,
                                                                        dbsnp_vcf_name) + '\n')
        else:
                print('exists')
                
        #Скачивание архивов dbSNP
        #из хранилища Ensembl.
        #До компьютера исследователя
        #они дойдут не целиком.
        #Пайплайну, получающему
        #короткие варианты, нужны
        #будут только хромосомы,
        #позиции, rsIDs и аллели.
        #По хромосоме и позиции
        #завершающая часть снипового
        #пайплайна будет искать rsID,
        #а по аллелям придётся иногда
        #доотбирать подходящую инсерцию.
        #Tabix-индексация поредевших таблиц.
        dbsnp_tsv_paths = []
        with open(dbsnp_urlstxt_path) as dbsnp_urlstxt_opened:
                for line in dbsnp_urlstxt_opened:
                        dbsnp_vcf_url = line.rstrip()
                        chr_name = re.search(r'(?<=chr)\w+(?=\.vcf\.gz)',
                                             dbsnp_vcf_url).group()
                        dbsnp_tsv_name = f'{chr_name}.tsv.gz'
                        dbsnp_tsv_path = os.path.join(ref_dir_path,
                                                      dbsnp_tsv_name)
                        dbsnp_tsv_paths.append(dbsnp_tsv_path)
                        
                        print(f'\n{dbsnp_tsv_name}', end='... ')
                        if os.path.exists(dbsnp_tsv_path) == False:
                                run_command(f'''
curl -s {dbsnp_vcf_url} |
zgrep "^[^#]" |
cut -f 1-5 |
bgzip -l 9 -@{threads_quan} > {dbsnp_tsv_path}''')
                        else:
                                print('exists')
                                
                        print(f'{dbsnp_tsv_name}.tbi', end='... ')
                        if os.path.exists(f'{dbsnp_tsv_path}.tbi') == False:
                                run_command(f'tabix -p vcf {dbsnp_tsv_path}')
                        else:
                                print('exists')
                                
        return dbsnp_tsv_paths
