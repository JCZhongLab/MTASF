For manual and reference, please visit http://proteomics.informatics.iupui.edu/software/toppic/

## System requirements

* GCC version higher than 4.8.2 for C++11 support
* CMake (>= 3.1)

### Windows:

[MSYS2](http://www.msys2.org/) is used for Windows building. Please follow the instructions from [here](doc/windows_build.md).

[MSYS2](http://www.msys2.org/) provides an excellent building platform for Windows. Please follow the instructions on MSYS2 website for installation and update the package database and core system packages. Then we can install other required packages using the commands below in MSYS2 shell:

```sh
pacman -S git

pacman -S mingw-w64-x86_64-gcc

pacman -S mingw-w64-x86_64-make

pacman -S mingw-w64-x86_64-cmake

pacman -S mingw-w64-x86_64-boost

pacman -S mingw-w64-x86_64-qt5

pacman -S mingw-w64-x86_64-xerces-c

pacman -S mingw-w64-x86_64-xalan-c
```

Thanks to MSYS2, we can install them very easily. Now we have C:\msys64\mingw64\include and C:\msys64\mingw64\lib for include and linking.

After installing, please add C:\msys64\mingw64\bin into your PATH environmental variable. In the following instructions, we will use both the MSYS2 shell and Windows CMD.

**Please make sure you can use** **g++** **in both MSYS2 shell and Windows CMD.**

**We use Windows CMD to build TopPIC suite.**

```sh
cd toppic-suite-MTASF
mkdir build
cd build
cmake -G "MinGW Makefiles" ..
mingw32-make
```

#### Attention:

If the file of xalan cannot be found during the compilation process, find the file name starting with libxalan, if you find the lib file whose name starts with liblibxalan, delete the extra string lib.

Modify the code using any IDE can be, it is recommended to use VS2019.

After compilation, copy the toppic_resources folder to the bin directory and compile the following instructions. Subsequent compilation only needs to execute mingw32-make in the build folder.

```sh
cd toppic-suite-MTASF
mkdir build
cd build
cmake -G "MinGW Makefiles" ..
mingw32-make
```

#### Example of input data:

- human_proteome_database.fasta
- RP4H_P32_WHIM2_biorep1_techrep1_ms2.msalign
- WHIM12_topmg_mods.txt

#### Test cases:

- ##### Open cmd or powershell in the bin directory and enter the following command:

```
topmg -x -i WHIM12_topmg_mods.txt -d -t FDR -v 0.01 -T FDR -V 0.01 human_proteome_database.fasta RP4H_P32_WHIM2_biorep1_techrep1_ms2.msalign
```

- ##### Enter the initial filtered database size  *f* , proteoform truncation ratio *k*, tag error tolerance range (ErrorT), and tag length *L* in turn according to the prompt.

  ![1650943575330](C:\Users\zyc\AppData\Roaming\Typora\typora-user-images\1650943575330.png)





### Linux (Ubuntu):

```sh
# install compiling tools
sudo apt-get install build-essential cmake

# install other dependencies
sudo apt-get install zlib1g-dev libboost-filesystem-dev \
                     libboost-program-options-dev \
                     libboost-system-dev \
                     libboost-thread-dev \
                     libboost-iostreams-dev \
                     libboost-chrono-dev \
                     libxalan-c-dev

# install the catch unit test framework (https://github.com/philsquared/Catch)
sudo apt-get install catch

# Qt5 for GUI
sudo apt-get install qtbase5-dev

# building
mkdir build
cd build
cmake ..
make -j$(nproc)

cd ../bin
ln -s ../toppic_resources .
```

On some Linux distributions, you might meet the problem "Could not loading a transcoding service".
To fix this, please add following lines into your `.bashrc`.

```sh
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
export LANGUAGE=en_US.UTF-8
```

### 