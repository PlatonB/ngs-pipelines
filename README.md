# Синопсис.
## Преимущества.
- Активно разрабатываемые компоненты от известных институтов и фирм.
- Высокая производительность:
  - Настраиваемое количество одновременно используемых ядер процессора;
  - По возможности задействуется Pipe для параллельного выполнения команд.
- Простота запуска.

## Перед началом работы.
### Установка сторонних компонентов.
Здесь придётся чуть покрасноглазить ☺.

#### Преодолеваем страх командной строки Linux.
Это лирическое отступление написано не для биоинформатиков, а для медиков/биологов, отважившихся обрабатывать геномные данные своими усилиями. [Пропустить↓](#установка-conda)

Основы основ:
- В Linux программы, позволяющие запускать команды, называются эмуляторами терминала. В составе вашего дистрибутива таковой будет, скорее всего, под названием **_Терминал_** или **_Консоль_**. Для тру-хакеров ☺ есть ещё чисто текстовое окружение _tty_, в которое можно переключиться сочетанием клавиш `CTRL+ALT+F2`. В народе любой инструмент выполнения команд называют терминалом, консолью и командной строкой.
- Процесс ввода и запуска команд ничем не отличаются от отправки сообщений в каком-нибудь мессенджере ☺. Просто наберите (или вставьте) команду и нажмите `↵`.
- Синтаксис, как правило, такой:
```
имя_программы собственно_команда --опция_1 --опция_2 --опция_N -i путь_к_исходному_файлу -o путь_к_конечному_файлу
```

В некоторых случаях, например, при запуске программы в _Docker_, синтаксис будет сложнее.

- `CTRL+ALT+C` - скопировать что-либо из эмулятора терминала (в [_elementary OS_](https://elementary.io/ru/) — проще — `CTRL+C`);
- `CTRL+ALT+V` - вставить что-либо в эмулятор терминала (в _elementary OS_ — проще — `CTRL+V`).
- `CTRL+D` - завершить работу запущенной в эмуляторе терминала программы. Это как нажатие на всем знакомый крестик у окон в графических окружениях.
- `sudo` перед командой нужен для того, чтобы совершить планируемое действие как суперпользователь (root). Обычно из-под рута проводятся рискованные манипуляции, требующие знаний и опыта. При работе на собственном компьютере вы можете легко задать root-пароль либо во время установки дистрибутива, либо уже потом командой `sudo passwd`. При вводе или вставке пароля вы, возможно, будете ожидать появление звёздочек. Но этого, вероятнее всего, не произойдёт: пароли в линуксовой командной строке, как правило, набираются/вставляются в невидимом режиме.
- `имя_программы --help` - вывести официальную документацию этой программы.
- `↑` - показать ранее выполненную команду. Ведь удобно же просмотреть (или заодно отредактировать и выполнить) одну из прошлых команд, не набирая её заново.

- Длинные команды удобно разбивать на несколько строк:
```
часть_команды \
другая_часть_команды \
и_так_далее
```
- `команда_1 && команда_2 && команда_N` - последовательное выполнение нескольких команд.
- `команда_1 | команда_2 | команда_N` - каждая команда передаёт вывод (stdout) следующей. Это позволяет одновременно (параллельно) выполнять задачи, что значительно ускоряет вычисления. Поддерживается не всеми консольными программами.

Часто пригождающиеся команды:
- `cd путь_к_папке` - перейти в нужную папку. Это позволит работать с файлами этой папки, не указывая абсолютный путь к каждому из них;
  - `cd ~` - перейти в папку home (домашняя папка пользователя);
  - `cd ..` - перейти в папку одним уровнем выше текущей.
- `ls` - вывести содержимое текущей папки.
- `mkdir путь_к_создаваемой_папке` - создать папку.
- `rm -r путь` - удалить папку со всеми вложенными объектами, либо файл.
- `имя_программы` - запустить программу, установленную с помощью пакетного менеджера вашего дистрибутива.
  - `./имя_программы` - запустить собранную вручную программу. Перед выполнением перейдите в папку, куда эта программа скомпилирована (`cd <...>`);
  - `sh имя_скрипта.sh` - выполнить скрипт, написанный на языке _Shell_. Не забывайте про `cd`;
  - `python3 имя_скрипта.py` - выполнить скрипт, написанный на языке Python 3. Тоже требуется `cd`.
- `less -N имя_файла_1.txt` - пролистывать файл колесом мыши. Если файл сжат _gzip_, то вместо `less` - `zless`. `-N` - нумерация строк. Выйти из просмотра - `q`.
- `grep 'поисковый_запрос' имя_файла_1.txt имя_файла_2.txt` - найти в файлах строки, содержащие запрашиваемое слово. Для поиска по gzip-архивам нужно применять `zgrep'.
- `htop` - запуск консольного монитора процессов (если тот установлен в ОС). При запуске громоздких задач (тех же NGS-пайплайнов) важно отслеживать использование оперативной памяти, нагрузку на ядра процессора и т.д..

Эмулятор терминала всё ещё кажется территорией суровых админов? ☺ Разъясню различные тонкости в [Issues](https://github.com/PlatonB/ngs-pipelines/issues).

#### Установка Conda.
1. Качаем [установочный скрипт](https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh).
2. `cd путь_к_папке_с_миникондой && sh Miniconda3-latest-Linux-x86_64.sh` - запускаем его.
3. Отвечаем `yes` на все задаваемые скриптом вопросы.
4. Подключаем нужные _Conda_-каналы:
```
conda config --add channels defaults
```
```
conda config --add channels bioconda
```
```
conda config --add channels conda-forge
```

#### Установка зависимостей.
Все сторонние компоненты, кроме _DeepVariant_, ставим через _Conda_. _Pysam_ потребуется только, если вы задумали обрабатывать человеческие данные.
```
conda install bowtie2 samtools pysam
```

_Google DeepVariant_ поставляется в виде _Docker_-образа, поэтому установка (на примере _Ubuntu_) будет такой:
```
sudo apt install docker.io
```
```
sudo docker pull google/deepvariant
```

Чтобы _Docker_ работал без навязчивых вводов root-пароля, надо выполнить несколко дополнительных действий.
```
sudo groupadd docker
```
```
sudo usermod -aG docker $USER
```
```
newgrp docker
```

В конце желательно перезагрузиться. Это, кстати, тоже можно из командной строки сделать.
```
systemctl reboot
```

Очень важно зависимости обновлять. От этого пайплайны будут работать быстрее и стабильнее. Запускайте команду обновления хотя бы раз в пару месяцев.
```
conda update --all
```

#### Запуск.
Запустить тот или иной пайплайн с помощью среды _IDLE_, установленной штатным менеджером пакетов, к сожалению, не получится. Но командная строка нам [теперь привычна↑](#преодолеваем-страх-командной-строки-linux), поэтому воспользуемся её услугами. Краткая схема построения команды:
```
python путь_к_пайплайну/имя_пайплайна.py обязательные_аргументы необязательные_аргументы
```

Все опции и соответствующие русскоязычные подсказки можно вывести так:
```
python путь_к_пайплайну/имя_пайплайна.py -h
```

#### Работа на сервере.
Мои программы успешно запускаются на обычном компе, но для их использования в реальной медицинской практике, конечно же, понадобится сервер.

В базовом случае, подключиться к серверу можно так:
```
ssh логин@адрес
```

Какой логин и куда подключаться, вам расскажет админ сервера.

Перекидывать со своего компа тексты программ и исходные данные можно с помощью обычных файлменеджеров: _Nautilus_, _Dolphin_, [_Pantheon Files_](https://elementary.io/ru/)... В последнем из перечисленных, к примеру, это делается нажатием `Подключить сервер` в боковом меню.

Для установки _Conda_ и удовлетворения зависимостей см. [инструкции выше↑](#установка-conda).

Если просто взять, запустить на сервере длительные вычисление и уйти надолго от компа, произойдёт разрыв канала связи. Выход из положения - применять утилиту _tmux_.

При небольшом количестве сэмплов и нежелании заморачиваться с _GNU Parallel_ создавайте для каждой задачи свою сессию. Чтобы не запутаться, лучше давать сессиям удобочитаемые имена. Допустим, в формате `имя_программы__фио_клиента__тип_секвенирования`.
```
tmux new -s variants_pipeline__platonb__wes
```

Войти в какую-либо сессию можно таким образом:
```
tmux attach -t variants_pipeline__platonb__wes
```

[Запустите↑](#запуск) внутри сессии нужную программу.

Пока идут подсчёты, ни в коем случае не выходите из сессии с помощью `CTRL+D`! Сессию надо аккуратно свернуть:
```
CTRL+B
```
и после этого
```
D
```

В свёрнутой сессии программа продолжит фурычить даже после разлогина. А разлогиниться можно уже обычным способом:
```
CTRL+D
```

Вновь залогинившись, можно в любой момент опять подсоединиться к сессии с помощью `tmux attach <...>`

Напоследок, пара полезных советов:
- `tmux ls` - вывести список сессий.
- `CTRL+B`, `[` - активировать режим прокрутки окна колесом мыши или двумя пальцами по тачпаду. Иначе вы будете пролистывать не вывод программы, а список прошлых команд. Выйти из этого режима - `esc`.

# variants_pipeline.
| Компонент | Выполняемая в рамках пайплайна работа |
| --------- | ------------------------------------- |
| _Bowtie 2_ | Индексация референсного генома, выравнивание |
| _BGZIP_	| Сжатие референсного генома |
| _SAMtools_ | Создание fai-индекса референсного генома, конвертация SAM в BAM, сортировка BAM, создание bai-индекса BAM |
| _DeepVariant_ | Коллинг вариантов |
| _Pysam_ | Дополнение VCF-файла идентификаторами вариантов. Только для Homo Sapiens |

Пример команды получения вариантов и их характеристик по непарным ридам полного человеческого экзома в 8 потоков.
```
python variants_pipeline.py \
-S $HOME/ngs/unpaired_reads \
-G $HOME/ngs/hg38_ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
-O homo_sapiens \
-R unpaired \
-t $HOME/ngs/trg \
-p 8 \
-s wes
```
